#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include <unistd.h> // read
#include <cstring>  // memchr
#include <algorithm>
#include <math.h>     
#include <omp.h> 
#include <sys/stat.h> // Check file exists

#include <tsl/robin_map.h>


#include <thread>

#include <chrono> // Timing functions

// g++ -std=c++11 -O3 -I$HOME/local/include PCG_BFA_4.cpp -o PCG_BFA_4 -fopenmp

using namespace std;

static double diffclock(clock_t clock1,clock_t clock2)
{
    double diffticks=clock1-clock2;
    double diffms=(diffticks)/(CLOCKS_PER_SEC/1000);
    return diffms;
}

int file_exists(char *name){
  struct stat   buffer;
  return (stat (name, &buffer) == 0);
}
int string_to_int(const char *p) {
    int x = 0;
    bool neg = false;
    if (*p == '-') {
        neg = true;
        ++p;
    }
    while (*p >= '0' && *p <= '9') {
        x = (x*10) + (*p - '0');
        ++p;
    }
    if (neg) {
        x = -x;
    }
    return x;
}

class StringRef
{
private:
    char const*     begin_;
    int             size_;

public:
    int size() const { return size_; }
    char const* begin() const { return begin_; }
    char const* end() const { return begin_ + size_; }

    StringRef( char const* const begin, int const size )
        : begin_( begin )
        , size_( size )
    {}
};

vector<StringRef> split4(const char *str, int size, char delimiter = '	'){
	vector<StringRef> results;
	const char *p = str;
	const char *q = str;
	
	while(p = (char*) memchr(p,delimiter,(str+size) - p)){
		results.emplace_back(StringRef(q,p-q));

		++p;
		q = p;
	}
	results.emplace_back(StringRef(q,(str+size)-q));

	return results;
}

double calculate_likelihood(tsl::robin_map<string,double> &selection_coefficients, vector<tsl::robin_map<string,double>> &initial_frequencies,vector<vector<int>> &timepoints,vector<tsl::robin_map<string,vector<int>>> &counts, vector<string> &all_lineage_IDs, int &number_of_replicates, tsl::robin_map<string,double> &conjugate_gradient_s_hash, vector< tsl::robin_map<string,double>> &conjugate_gradient_f0_hash, double &alpha){
	long double log_likelihood = 0;

	int number_of_lineages = all_lineage_IDs.size();
	for(int rep = 0;rep < number_of_replicates;++rep){
		// Expected f(t) are replicate dependent, easier to loop this way.
		int timepoints_size = timepoints[rep].size();
		long double expected;
		
		for(int t = 0;t < timepoints_size;++t){
			long double total = 0;
			long double sum = 0;
			for(int i = 0;i < number_of_lineages;++i){
				if(counts[rep].find(all_lineage_IDs[i]) != counts[rep].end()){
					// Expected f(t) = f(0) * exp(st)
					long double f0 = initial_frequencies[rep][all_lineage_IDs[i]] + alpha * conjugate_gradient_f0_hash[rep][all_lineage_IDs[i]];
					long double s = selection_coefficients[all_lineage_IDs[i]] + alpha * conjugate_gradient_s_hash[all_lineage_IDs[i]];
					
					expected = f0 * exp(s * (long double) timepoints[rep][t]) ;
					sum  += expected ;
					
					log_likelihood += (long double) counts[rep][all_lineage_IDs[i]][t] * log(expected);
					total += (long double) counts[rep][all_lineage_IDs[i]][t];
				}
			}
			log_likelihood -= (long double) total * log(sum);

		}
		
	}
	return log_likelihood;
}
double derivative_wrt_alpha(tsl::robin_map<string,double> &selection_coefficients, vector<tsl::robin_map<string,double>> &initial_frequencies,vector<vector<int>> &timepoints,vector<tsl::robin_map<string,vector<int>>> &counts, vector<vector<int>> &total_reads, vector<string> &all_lineage_IDs, int &number_of_replicates, tsl::robin_map<string,double> &conjugate_gradient_s_hash, vector< tsl::robin_map<string,double>> &conjugate_gradient_f0_hash, double &alpha){
	long double derivative = 0;

	int number_of_lineages = all_lineage_IDs.size();
	for(int rep = 0;rep < number_of_replicates;++rep){
		// Expected f(t) are replicate dependent, easier to loop this way.
		int timepoints_size = timepoints[rep].size();
		
		for(int t = 0;t < timepoints_size;++t){
			long double total = 0;
			long double sum_top = 0;
			long double sum_bottom = 0;
			
			for(int i = 0;i < number_of_lineages;++i){
				if(counts[rep].find(all_lineage_IDs[i]) != counts[rep].end()){
					int nk = counts[rep][all_lineage_IDs[i]][t];
					long double ck = conjugate_gradient_f0_hash[rep][all_lineage_IDs[i]];
					long double phik = conjugate_gradient_s_hash[all_lineage_IDs[i]];
					// Expected f(t) = f(0) * exp(st)		
					long double f0 = initial_frequencies[rep][all_lineage_IDs[i]] + alpha * conjugate_gradient_f0_hash[rep][all_lineage_IDs[i]];
					long double s = selection_coefficients[all_lineage_IDs[i]] + alpha * conjugate_gradient_s_hash[all_lineage_IDs[i]];
					
					long double est =  exp(f0+ s * timepoints[rep][t] );
					sum_bottom += est;
					sum_top += total_reads[rep][t] * est * (ck + phik * timepoints[rep][t]);
					
					derivative += (ck + phik * timepoints[rep][t]) * counts[rep][all_lineage_IDs[i]][t];
				}
			}
			
			derivative -= sum_top/sum_bottom;
			
		}
	}
	//cout << alpha << "	" << derivative << endl;
	
	return derivative;
}
std::pair<double, double> alpha_newton_step(tsl::robin_map<string,double> &selection_coefficients, vector<tsl::robin_map<string,double>> &initial_frequencies,vector<vector<int>> &timepoints,vector<tsl::robin_map<string,vector<int>>> &counts, vector<vector<int>> &total_reads, vector<string> &all_lineage_IDs, int &number_of_replicates, tsl::robin_map<string,double> &conjugate_gradient_s_hash, vector< tsl::robin_map<string,double>> &conjugate_gradient_f0_hash, double &alpha){
	double derivative = 0;
	double second_derivative = 0;
	int number_of_lineages = all_lineage_IDs.size();

	for(int rep = 0;rep < number_of_replicates;++rep){
		// Expected f(t) are replicate dependent, easier to loop this way.
		int timepoints_size = timepoints[rep].size();
		
		#pragma omp parallel for reduction(+:second_derivative, derivative)
		for(int t = 0;t < timepoints_size;++t){
			double sum_top_left = 0;
	
			double sum_top = 0;
			double sum_bottom = 0;

			// Need to prevent overflow of the exponential function. We'll try the exp normalize trick.
			// First let's obtain the maximum value. There must be a way by omp to get this.
			double max_val = 0;
			#pragma omp parallel for reduction(max: max_val)
			for(int i = 0;i < number_of_lineages;++i){
				if(counts[rep].find(all_lineage_IDs[i]) != counts[rep].end()){
					double f0 = initial_frequencies[rep][all_lineage_IDs[i]] + alpha * conjugate_gradient_f0_hash[rep][all_lineage_IDs[i]];
					double s = selection_coefficients[all_lineage_IDs[i]] + alpha * conjugate_gradient_s_hash[all_lineage_IDs[i]];

					if(f0 + s * timepoints[rep][t] > max_val){
						max_val = f0 + s * timepoints[rep][t];
					}
				}
			}


			#pragma omp parallel for reduction(+:sum_bottom, sum_top, sum_top_left)
			for(int i = 0;i < number_of_lineages;++i){
				if(counts[rep].find(all_lineage_IDs[i]) != counts[rep].end()){
					//printf("thread is %d\n", omp_get_thread_num());
					int nk = counts[rep][all_lineage_IDs[i]][t];
					//cout << rep << "	" << i << "	" << all_lineage_IDs[i] << "	" << conjugate_gradient_f0_hash[rep][all_lineage_IDs[i]] << endl;
					double ck = conjugate_gradient_f0_hash[rep][all_lineage_IDs[i]];
					double phik = conjugate_gradient_s_hash[all_lineage_IDs[i]];
					
					// Expected f(t) = f(0) * exp(st)		
					double f0 = initial_frequencies[rep][all_lineage_IDs[i]] + alpha * conjugate_gradient_f0_hash[rep][all_lineage_IDs[i]];
					double s = selection_coefficients[all_lineage_IDs[i]] + alpha * conjugate_gradient_s_hash[all_lineage_IDs[i]];
					
					double est =  exp(f0+ s * timepoints[rep][t] - max_val);
					double ck_phik_t = (ck + phik * timepoints[rep][t]);
					
					sum_bottom += est;
					sum_top += est * ck_phik_t;
					
					derivative += ck_phik_t * counts[rep][all_lineage_IDs[i]][t];	
					sum_top_left += est * ck_phik_t * ck_phik_t;
				}
			}
			derivative -= total_reads[rep][t] * sum_top/sum_bottom;
			second_derivative += (total_reads[rep][t] * (sum_top_left / sum_bottom - sum_top * sum_top / (sum_bottom*sum_bottom)));
		}
	}
	second_derivative *= -1;
	
	return std::make_pair(derivative,second_derivative);
}
int main(int argc, char * argv[]) {
	cin.tie(NULL);

	std::cout.unsetf ( std::ios::floatfield ); 
	std::cerr.unsetf (std::ios::floatfield );  
	std::cerr.precision(20);
	std::cout.precision(20);

	omp_set_num_threads(omp_get_num_procs());

	//cout << omp_get_num_procs() << endl;
	istream *in;
	// Must be declared here for scope reasons
	ifstream ifn;
	
	// Read the frequencies files. If there are multiple files, read all of them.
	int buffer_size = 512;
	char * buffer = new char [buffer_size];
	
	int number_of_replicates = argc-1;
	
	// Generate the massive array of data.
	vector<vector<int>> timepoints;
	vector<vector<int>> total_reads;
	vector<tsl::robin_map<string,vector<int>>> counts;
	vector<string> all_lineage_IDs;
	
	tsl::robin_map<string,int> tmp_lineage_IDs;
	tsl::robin_map<string,int> sum_lineage_IDs;
	
	tsl::robin_map<string,double> selection_coefficients;
	vector<tsl::robin_map<string,double>> initial_frequencies;
	
	if(number_of_replicates < 1){
		cerr << "gradient_BFA_reps file [file] [file]" << endl;
		return 1;
	}

	for(int i = 1;i < argc;++i){
		if(file_exists(argv[i])){
			string prev_line = "";
				
			timepoints.push_back(vector<int>());
			total_reads.push_back(vector<int>());
			counts.push_back(tsl::robin_map<string,vector<int>>());
			initial_frequencies.push_back(tsl::robin_map<string,double>());
			
			ifn.open(argv[i]);
			in = &ifn;
			
			for(string row; getline( *in, row); ){
				vector<StringRef> const fields = split4(row.c_str(),row.length());
					
				// fields now represent a vector of the tab delimited file.
				// We now parse the first row, since it contains the header information about the timepoints.
				if(timepoints[i-1].size() == 0){
					for(int j = 1;j < fields.size();++j){
						timepoints[i-1].push_back(string_to_int(string(fields[j].begin(),fields[j].end()).c_str()));
						total_reads[i-1].push_back(0);
					}
				}
				else{
					// Now we have already parsed the first row, so we want to store the read counts.
					// The lineage ID is the first column.
					// First check if it's present (all counts > 0)
					int present = 0;
					for(int j = 1;j < fields.size();++j){
						if(string_to_int(string(fields[j].begin(),fields[j].end()).c_str()) > 0){
							present = 1;
							break;
						}
					}
						
					if(present == 1){
						string lineage_ID = string(fields[0].begin(),fields[0].end()).c_str();
						// Store the lineage ID
						tmp_lineage_IDs[lineage_ID] += 1;
						
						for(int j = 1;j < fields.size();++j){
							// Grab the counts, and add up the counts to the total
							int this_count = string_to_int(string(fields[j].begin(),fields[j].end()).c_str());							
							counts[i-1][lineage_ID].push_back(this_count);
							sum_lineage_IDs[lineage_ID] += this_count;
							
							total_reads[i-1][j-1] += this_count;
						}
					}
				}

			}

			ifn.close();
			ifn.clear();
		}
		else{
			cerr << "Cannot access file " << argv[i] << endl;
			return 1;
		}
	}
	

	// Done reading files
	// Now find the reference lineage.
	// To do so, we will use a lineage that is found in all the replicates, with the total highest counts.
	// Let's fill the array of all the lineage IDs.
	int max_count = 0;
	int iter_pos_ref_lineage = -1;
	for(auto it = tmp_lineage_IDs.begin(); it != tmp_lineage_IDs.end(); ++it) {
		all_lineage_IDs.push_back(it->first);
		// Ask if it's found in all the replicates.
		if(it->second == argc-1 && sum_lineage_IDs[it->first] > max_count){
			max_count = sum_lineage_IDs[it->first];
			
			iter_pos_ref_lineage = all_lineage_IDs.size();
		}
	}

	if(iter_pos_ref_lineage == -1){
		cerr << "Could not find a reference lineage." << endl;
		return 1;
	}

	int number_of_lineages = all_lineage_IDs.size();


	// Reference lineage currently at: all_lineage_IDs[iter_pos_ref_lineage-1]
	// Move to the end by swapping position with the last element. 
	iter_swap(all_lineage_IDs.end()-1, all_lineage_IDs.begin() + iter_pos_ref_lineage-1);

	// Now we have completely parsed the files.
	// We must now initialize the regression values. These are the parameters we are estimating.
	for(int i = 0;i < number_of_lineages;++i){
		selection_coefficients[all_lineage_IDs[i]] = 0;
		
		#pragma omp parallel for		
		for(int rep = 0;rep < number_of_replicates;++rep){
			// First we have to make sure that the lineage exists in this replicate.
			if(counts[rep].find(all_lineage_IDs[i]) != counts[rep].end()){
				initial_frequencies[rep][all_lineage_IDs[i]] = (1.0+counts[rep][all_lineage_IDs[i]][0])/(counts[rep].size() + total_reads[rep][0]);
			}
		}
	}

	// Done initializing values.
	// Now obtain a reasonable guess.
	double prev_log_likelihood = -1*std::numeric_limits<double>::infinity();;
	for(int iteration = 0;iteration < 100;++iteration){
		// We do as few iterations since multiple iterations aren't guaranteed to increase the likelihood.
		double log_likelihood = 0.0;
		// First, get the likelihood of the 'naive' initialization.
		// The likelihood will be calculated as a multinomial. We will ignore the constant, which is not necessary for most purposes.
		// sum n * log (f (t))
		
		// We must first get the expected.
		vector< vector<double>> sum_expected;
		vector< tsl::robin_map<string,vector<double>>> expected;
		sum_expected.resize(number_of_replicates);
		for(int rep = 0;rep < number_of_replicates;++rep){
			// Expected f(t) are replicate dependent, easier to loop this way.
			sum_expected[rep].resize(timepoints[rep].size(),0);
			expected.push_back(tsl::robin_map<string,vector<double>>());
			
			for(int t = 0;t < timepoints[rep].size();++t){
				for(int i = 0;i < number_of_lineages;++i){
					if(counts[rep].find(all_lineage_IDs[i]) != counts[rep].end()){
						// Expected f(t) = f(0) * exp(st)
						expected[rep][all_lineage_IDs[i]].push_back( initial_frequencies[rep][all_lineage_IDs[i]] * exp(selection_coefficients[all_lineage_IDs[i]] * timepoints[rep][t]) );
						sum_expected[rep][t] += expected[rep][all_lineage_IDs[i]][t];
					}
				}
				
				// We now have the sum expected, so we normalize.
				//#pragma omp parallel for reduction(+:log_likelihood)
				for(int i = 0;i < number_of_lineages;++i){
					if(counts[rep].find(all_lineage_IDs[i]) != counts[rep].end()){
						// Expected f(t) = f(0) * exp(st)	/ sum
						expected[rep][all_lineage_IDs[i]][t] /= sum_expected[rep][t];

						// Now calculate the log likelihood
						log_likelihood += counts[rep][all_lineage_IDs[i]][t] * log(expected[rep][all_lineage_IDs[i]][t]);
					}
				}
			}
		}
		
		cerr << "Log likelihood at the beginning of iteration: " << log_likelihood << endl;
		if(log_likelihood <= prev_log_likelihood){
			// Have done our best with the standard regression. Now we'll optimize finely.
			break;
		}
		prev_log_likelihood = log_likelihood;
		
		vector< tsl::robin_map<string,double>> intercepts;
		vector<double> sum_intercepts;
		intercepts.resize(number_of_replicates);
		sum_intercepts.resize(number_of_replicates);
		
		// Done calculating likelihood, now perform log linear regression.
		for(int i = 0;i < number_of_lineages-1;++i){
		
			vector<double> sum_w;
			vector<double> sum_xx;
			vector<double> sum_xy;
			vector<double> sum_x;
			vector<double> sum_y;
			
			sum_w.resize(number_of_replicates);
			sum_xx.resize(number_of_replicates);
			sum_xy.resize(number_of_replicates);
			sum_x.resize(number_of_replicates);
			sum_y.resize(number_of_replicates);
			
			double slope_up = 0;
			double slope_down = 0;
		
			for(int rep = 0;rep < number_of_replicates;++rep){
				if(counts[rep].find(all_lineage_IDs[i]) != counts[rep].end()){
					
					for(int t = 0;t < timepoints[rep].size();++t){
						double expected_val = expected[rep][all_lineage_IDs[i]][t];
						double expected_zero = expected[rep][all_lineage_IDs[number_of_lineages-1]][t];
						
						// Obtain the expected variance.
						double x_1 = sqrt((1-expected_val) / (total_reads[rep][t] * expected_val));
						double x_2 = sqrt((1-expected_zero) / (total_reads[rep][t] * expected_zero));
						
						double var_1 = x_1 - x_1*x_1/2;
						double var_2 = x_2 - x_2*x_2/2;
						
						var_1 = var_1*var_1;
						var_2 = var_2*var_2;

						double variance = var_1+var_2;
						
						double w = 1/(variance);
						
						sum_w[rep] += w;
						
						// Now we obtain the observed fraction so that we may regress on it.
						double count_F = (double) counts[rep][all_lineage_IDs[i]][t]/total_reads[rep][t];
						double count_Z = (double) counts[rep][all_lineage_IDs[number_of_lineages-1]][t]/total_reads[rep][t];
						
						if(count_F == 0){
							count_F = (1.0-pow(0.5,1.0/total_reads[rep][t]))/2.0;
						}
						if(count_Z == 0){
							count_Z = (1.0-pow(0.5,1.0/total_reads[rep][t]))/2.0;
						}
						
						
						
						double y = log(count_F/count_Z);
						double x = timepoints[rep][t];
					
						
						sum_xy[rep] += w*x*y;
						sum_xx[rep] += w*x*x;
						sum_x[rep] += w*x;
						sum_y[rep] += w*y;
						
					}
					
					slope_up += sum_xy[rep] - sum_y[rep]*sum_x[rep]/sum_w[rep];
					slope_down += sum_xx[rep] - sum_x[rep]*sum_x[rep]/sum_w[rep]; 
					
				}
			}
			
			// Now calculate the selection coefficient.
			selection_coefficients[all_lineage_IDs[i]] = slope_up/slope_down;
					
			// Calculate the initial frequencies now.
			for(int rep = 0;rep < number_of_replicates;++rep){
				if(counts[rep].find(all_lineage_IDs[i]) != counts[rep].end()){
					double intercept = (sum_y[rep] - selection_coefficients[all_lineage_IDs[i]] * sum_x[rep])/sum_w[rep];

					intercept = exp(intercept);
					
					intercepts[rep][all_lineage_IDs[i]] = intercept;
					sum_intercepts[rep] += intercept;
				}
			}
		}
		
		// Now we make sure all intercepts sum to 1.
		for(int rep = 0;rep < number_of_replicates;++rep){
			double sum_not_ref = 0;
			
			#pragma omp parallel for reduction(+:sum_not_ref)
			for(int i = 0;i < number_of_lineages-1;++i){
				if(counts[rep].find(all_lineage_IDs[i]) != counts[rep].end()){
					initial_frequencies[rep][all_lineage_IDs[i]] = intercepts[rep][all_lineage_IDs[i]]/(1+sum_intercepts[rep]);	

					sum_not_ref += initial_frequencies[rep][all_lineage_IDs[i]];
				}
			}
			
			initial_frequencies[rep][all_lineage_IDs[number_of_lineages-1]] = 1.0 - sum_not_ref;
		}
	}
	
	// We now reparametrize the problem into the following:
	// f(0) = log(f(0)/fref(0))
	// Such that fref = 0
	for(int rep = 0;rep < number_of_replicates;++rep){
		#pragma omp parallel for
		for(int i = 0;i < number_of_lineages-1;++i){
			if(counts[rep].find(all_lineage_IDs[i]) != counts[rep].end()){
				initial_frequencies[rep][all_lineage_IDs[i]] = 	log(initial_frequencies[rep][all_lineage_IDs[i]]/initial_frequencies[rep][all_lineage_IDs[number_of_lineages-1]]);

			}
		}
		
		initial_frequencies[rep][all_lineage_IDs[number_of_lineages-1]] = 0;
	}
	
	
	// Now we have our initial guess. Optimize further by conjugate gradient.
	double golden = 1.0-1.0/((1.0+sqrt(5.0))/2.0);
	tsl::robin_map<string,double> conjugate_gradient_s_hash;  // Replicate independent
	vector< tsl::robin_map<string,double>> conjugate_gradient_f0_hash;  // Replicate dependent
	conjugate_gradient_f0_hash.resize(number_of_replicates);
	
	tsl::robin_map<string,double> previous_gradient_s_hash;  // Replicate independent
	vector< tsl::robin_map<string,double>> previous_gradient_f0_hash;  // Replicate dependent
	previous_gradient_f0_hash.resize(number_of_replicates);
	double log_likelihood;
	
	double time_0, time_1;
	chrono::high_resolution_clock::time_point walltime_0, walltime_1;
	time_0 = clock();
	walltime_0 = chrono::high_resolution_clock::now();
	double tmp_value;
	int iteration_count;

	int reset_flag = 0;

	for(int iteration = 0;iteration < 100000;++iteration){
		
		log_likelihood = 0.0;

		// We must first get the expected.
		// This takes 5% of the computational time.
		vector< vector<double>> sum_expected;
		vector< tsl::robin_map<string,vector<double>>> expected;
		
		sum_expected.resize(number_of_replicates);
		expected.resize(number_of_replicates);
			
		for(int rep = 0;rep < number_of_replicates;++rep){
			int timepoints_size = timepoints[rep].size();
			sum_expected[rep].resize(timepoints_size,0);	
			for(int t = 0;t < timepoints_size;++t){	

				for(auto it = counts[rep].begin();it != counts[rep].end();++it){
					tmp_value = initial_frequencies[rep][it->first] + selection_coefficients[it->first] * timepoints[rep][t];
					expected[rep][it->first].push_back( exp(tmp_value) );
					sum_expected[rep][t] += expected[rep][it->first][t];
					
					log_likelihood += counts[rep][it->first][t] * tmp_value;
				}

				

				#pragma omp parallel for
				for(int i = 0;i < number_of_lineages-1;++i){
					if(counts[rep].find(all_lineage_IDs[i]) != counts[rep].end()){
						expected[rep][all_lineage_IDs[i]][t] /= sum_expected[rep][t];	
					}
				}
				
				log_likelihood -= total_reads[rep][t] * log(sum_expected[rep][t]);
			}
		}


		
		
		// Obtain the gradient and the conjugate gradient.
		// 3% of the computational time.
		tsl::robin_map<string,double> gradient_s_hash;  // Replicate independent
		vector< tsl::robin_map<string,double>> gradient_f0_hash;  // Replicate dependent
		gradient_f0_hash.resize(number_of_replicates);

		tsl::robin_map<string,double> hessian_s_hash;  // Replicate independent
		vector< tsl::robin_map<string,double>> hessian_f0_hash;  // Replicate dependent
		hessian_f0_hash.resize(number_of_replicates);
		
		double beta = 0;

		double gTg = 0;
		double dTy = 0; // dtg - dtg_old
		double dTg = 0;
		
		double yTy = 0; // gTg - 2*gTg_old + g_oldTg_old
		double yTg = 0; // gTg - gTg_old
		double dTg_old = 0;
		double dTd = 0;

		// g = gradient
		// d = conjugate gradient
		// y = difference in gradient
		
		// To rewrite for openmp
		// Still totally unclear why a reduction on the 7 variables does not work. Probably because inserting into the unordered maps is concurrent and not allowed.
		// To make this work for openMP, I need to remove writing into the hash... can't see how i would do that.
		double max_g = 0;
		for(int i = 0;i < number_of_lineages-1;++i){
			double gradient_s = 0;
			double gradient_s2 = 0;
			for(int rep = 0;rep < number_of_replicates;++rep){
				if(counts[rep].find(all_lineage_IDs[i]) != counts[rep].end()){
					double gradient_f0 = 0;	
					double gradient_f02 = 0;
					int timepoints_size = timepoints[rep].size();
					for(int t = 0;t < timepoints_size;++t){	
						double tmp_value = (counts[rep][all_lineage_IDs[i]][t] - total_reads[rep][t] * expected[rep][all_lineage_IDs[i]][t]);
						double tmp_value2 = total_reads[rep][t] * expected[rep][all_lineage_IDs[i]][t] * (1-expected[rep][all_lineage_IDs[i]][t]);
						gradient_s += timepoints[rep][t] * tmp_value;						
						gradient_f0 += tmp_value;
						gradient_s2 += timepoints[rep][t] * timepoints[rep][t] * tmp_value2;
						gradient_f02 += tmp_value2;
						
					
					}					
					gradient_f0_hash[rep][all_lineage_IDs[i]] = gradient_f0;
					hessian_f0_hash[rep][all_lineage_IDs[i]] = -1*gradient_f02;
					if(fabs(gradient_f0) > max_g){max_g = fabs(gradient_f0);}
					
					gTg += gradient_f0*gradient_f0;
					dTy += conjugate_gradient_f0_hash[rep][all_lineage_IDs[i]] * (gradient_f0-previous_gradient_f0_hash[rep][all_lineage_IDs[i]]);
					dTg += conjugate_gradient_f0_hash[rep][all_lineage_IDs[i]] * gradient_f0;
					yTy += (gradient_f0-previous_gradient_f0_hash[rep][all_lineage_IDs[i]])* (-1/gradient_f02) * (gradient_f0-previous_gradient_f0_hash[rep][all_lineage_IDs[i]]);
					dTg_old += conjugate_gradient_f0_hash[rep][all_lineage_IDs[i]] * previous_gradient_f0_hash[rep][all_lineage_IDs[i]];
					dTd += conjugate_gradient_f0_hash[rep][all_lineage_IDs[i]] * conjugate_gradient_f0_hash[rep][all_lineage_IDs[i]] * gradient_f02;
					yTg += (gradient_f0-previous_gradient_f0_hash[rep][all_lineage_IDs[i]]) * (-1/gradient_f02) * gradient_f0;
						
				}	
			}
	
			gradient_s_hash[all_lineage_IDs[i]] = gradient_s;
			hessian_s_hash[all_lineage_IDs[i]] = -1*gradient_s2;
			
			gTg += gradient_s*gradient_s;
			dTy += conjugate_gradient_s_hash[all_lineage_IDs[i]] * (gradient_s-previous_gradient_s_hash[all_lineage_IDs[i]]);		
			dTg += conjugate_gradient_s_hash[all_lineage_IDs[i]] * gradient_s;		
			yTy += (gradient_s-previous_gradient_s_hash[all_lineage_IDs[i]])* (-1/gradient_s2) * (gradient_s-previous_gradient_s_hash[all_lineage_IDs[i]]);
			dTg_old += conjugate_gradient_s_hash[all_lineage_IDs[i]] * previous_gradient_s_hash[all_lineage_IDs[i]];
			dTd += conjugate_gradient_s_hash[all_lineage_IDs[i]] * conjugate_gradient_s_hash[all_lineage_IDs[i]] * gradient_s2;
			yTg += (gradient_s-previous_gradient_s_hash[all_lineage_IDs[i]]) * (-1/gradient_s2) * gradient_s;
			if(fabs(gradient_s) > max_g){max_g = fabs(gradient_s);}				
		}
	

		previous_gradient_s_hash = gradient_s_hash;
		previous_gradient_f0_hash = gradient_f0_hash;

		if(iteration == 0 || reset_flag == 1){
			// Conjugate gradient is the gradient initially.
			conjugate_gradient_s_hash = gradient_s_hash;
			conjugate_gradient_f0_hash = gradient_f0_hash;
			reset_flag = 0;
		}
		else{
			double trunc;
			beta = -1.0*(yTg/dTy - (yTy/dTy)*(dTg/dTy)); // Hager-Zhang conjugate gradient.

			//beta = -1.0 * (yTg/dTy - (3*yTy/dTy - dTy/dTd) * (dTg/dTy)); 
			//cout << beta << "	" << trunc << endl;
			trunc = -0.4*dTg_old/dTd; // As presented in Dai Kou (they use 0.5 here). Note that beta is almost always used.
			beta = (beta > trunc ? beta : trunc);
			
		}

		if(iteration % 50 == 0){
			time_1 = clock();
			walltime_1 = chrono::high_resolution_clock::now();
			chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(walltime_1 - walltime_0);
			cerr << "Log likelihood at the beginning of iteration " << iteration << ": " << log_likelihood << "	" << "Objective: " << sqrt(gTg) << "	" << "Max gradient: " << max_g << "	Clock: " << diffclock(time_1,time_0)/1000 <<  "	" << "Wallclock: " << time_span.count() << endl;
		}
		
		// Beta is now measured. We now find the step size.
		// Get the conjugate gradient.

		#pragma omp parallel for
		for(int i = 0;i < number_of_lineages-1;++i){
			conjugate_gradient_s_hash[all_lineage_IDs[i]] = (1/hessian_s_hash[all_lineage_IDs[i]]) * gradient_s_hash[all_lineage_IDs[i]] + beta * conjugate_gradient_s_hash[all_lineage_IDs[i]];
			
			for(int rep = 0;rep < number_of_replicates;++rep){
				if(counts[rep].find(all_lineage_IDs[i]) != counts[rep].end()){
					conjugate_gradient_f0_hash[rep][all_lineage_IDs[i]] = (1/hessian_f0_hash[rep][all_lineage_IDs[i]]) * gradient_f0_hash[rep][all_lineage_IDs[i]] + beta * conjugate_gradient_f0_hash[rep][all_lineage_IDs[i]];
				}
			}
		}

		// Deal with the reference lineage for openMP
		for(int rep = 0;rep < number_of_replicates;++rep){
			conjugate_gradient_f0_hash[rep][all_lineage_IDs[number_of_lineages-1]] = 0;
		}	
		conjugate_gradient_s_hash[all_lineage_IDs[number_of_lineages-1]] = 0;

		// Find step size. We're going to do a safeguarded Newton-Raphson. 
		// Bracket the root.

		double min_alpha = 0;
		double optimal_alpha = min_alpha;

		// Find the maximum. We'll do a series of newton steps until the derivative is positive.
		
		double max_alpha = -1;

		//cout << gradient_s_hash[all_lineage_IDs[number_of_lineages-1]] << "	" << conjugate_gradient_s_hash[all_lineage_IDs[number_of_lineages-1]] << endl;
		//cout << gradient_f0_hash[1][all_lineage_IDs[number_of_lineages-1]] << "	" << conjugate_gradient_f0_hash[1][all_lineage_IDs[number_of_lineages-1]] << endl;


		optimal_alpha = max_alpha;
		while(1){
			pair<double, double> alpha_step = alpha_newton_step(selection_coefficients, initial_frequencies, timepoints, counts, total_reads, all_lineage_IDs, number_of_replicates, conjugate_gradient_s_hash, conjugate_gradient_f0_hash, optimal_alpha);	
			//return 0;
			

			if(alpha_step.first > 0){
				max_alpha = optimal_alpha;
				break;
			}
			min_alpha = optimal_alpha;
			optimal_alpha *= 2;
		}
		
		optimal_alpha = (min_alpha + max_alpha)/2;

		// We now have the brackets (min to max). We'll do a newton raphson step or bisect if the step is outside of bound.
		int k = 0;
		while(k < 20){
			pair<double, double> alpha_step = alpha_newton_step(selection_coefficients, initial_frequencies, timepoints, counts, total_reads, all_lineage_IDs, number_of_replicates, conjugate_gradient_s_hash, conjugate_gradient_f0_hash, optimal_alpha);
			double newton_step = alpha_step.first/alpha_step.second;
			double putative_optimal_alpha = optimal_alpha - newton_step;

			// Decrease the bounds
			if(alpha_step.first < 0){
				min_alpha = optimal_alpha;
			}
			else{
				max_alpha = optimal_alpha;
			}

			
			if(putative_optimal_alpha > min_alpha || putative_optimal_alpha < max_alpha){
			//if((putative_alpha - min_alpha) * (putative_alpha - max_alpha) > 0){ 
				// out of bounds.
				// Bisect.
				optimal_alpha = (min_alpha + max_alpha)/2;
			}
			else{
				optimal_alpha -= newton_step;
			}
			++k;
			if(fabs(newton_step) < 1e-8 || fabs(max_alpha - min_alpha) < 1e-8){ // Machine epsilon is 1e-16
				break;
			}
		}

		if(optimal_alpha > 0 || isnan(optimal_alpha)){
			//reset_flag = 1;
			optimal_alpha = -0.1;
		}
		//if(optimal_alpha < -1){
		//	optimal_alpha = -1;
		//}

		//cout << optimal_alpha << endl;
		
		// Now apply the gradient.
		//double max_change = -1.0*std::numeric_limits<double>::infinity();
		#pragma omp parallel for
		for(int i = 0;i < number_of_lineages-1;++i){
			selection_coefficients[all_lineage_IDs[i]] += optimal_alpha * conjugate_gradient_s_hash[all_lineage_IDs[i]];
			for(int rep = 0;rep < number_of_replicates;++rep){
				if(counts[rep].find(all_lineage_IDs[i]) != counts[rep].end()){
					initial_frequencies[rep][all_lineage_IDs[i]] += optimal_alpha * conjugate_gradient_f0_hash[rep][all_lineage_IDs[i]];
						
				}
			}
			
		}
		
		//return 1;
		
		if(sqrt(gTg) < 1){
			// Convergence criteria
			iteration_count = iteration;

			time_1 = clock();
			walltime_1 = chrono::high_resolution_clock::now();
			chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(walltime_1 - walltime_0);
			cerr << "Log likelihood at the beginning of iteration " << iteration << ": " << log_likelihood << "	" << "Objective: " << sqrt(gTg) << "	" << "Max gradient: " << max_g << "	Clock: " << diffclock(time_1,time_0)/1000 <<  "	" << "Wallclock: " << time_span.count() << endl;
		


			break;
		}
		
	}
	
	
	// Done, now compute the Fisher information matrix.
	// Also compute the squared residuals.	
	log_likelihood = 0.0;
	double full_log_likelihood = 0;
	double chi_square = 0;
	int df = 0; // Full parameter space calculation
	
	// We must first get the expected.
	vector< vector<double>> sum_expected;
	vector< tsl::robin_map<string,vector<double>>> expected;
	
	sum_expected.resize(number_of_replicates);
	expected.resize(number_of_replicates);	
	for(int rep = 0;rep < number_of_replicates;++rep){
		int timepoints_size = timepoints[rep].size();
		sum_expected[rep].resize(timepoints_size,0);	
		for(int t = 0;t < timepoints_size;++t){
			for(int i = 0;i < number_of_lineages;++i){
				if(counts[rep].find(all_lineage_IDs[i]) != counts[rep].end()){
					
					expected[rep][all_lineage_IDs[i]].push_back( exp(initial_frequencies[rep][all_lineage_IDs[i]] + selection_coefficients[all_lineage_IDs[i]] * timepoints[rep][t]) );
					sum_expected[rep][t] += expected[rep][all_lineage_IDs[i]][t];
					
					log_likelihood += (double)counts[rep][all_lineage_IDs[i]][t] * log(expected[rep][all_lineage_IDs[i]][t]);
				}
			}

			for(int i = 0;i < number_of_lineages;++i){
				if(counts[rep].find(all_lineage_IDs[i]) != counts[rep].end()){
					expected[rep][all_lineage_IDs[i]][t] /= sum_expected[rep][t];	
					chi_square += (double)(counts[rep][all_lineage_IDs[i]][t]-expected[rep][all_lineage_IDs[i]][t]*total_reads[rep][t])*(counts[rep][all_lineage_IDs[i]][t]-expected[rep][all_lineage_IDs[i]][t]*total_reads[rep][t])/(expected[rep][all_lineage_IDs[i]][t]*total_reads[rep][t]);
					
					if(counts[rep][all_lineage_IDs[i]][t] > 0){
						full_log_likelihood += (double)counts[rep][all_lineage_IDs[i]][t] * log((double) counts[rep][all_lineage_IDs[i]][t]/(double)total_reads[rep][t]);
					}
					++df;
				}
			}
			log_likelihood -= (double)total_reads[rep][t] * log(sum_expected[rep][t]);
			--df;
		}
	}
	// The number of parameters in our model is 2 * (number_of_lineages - 1).
	// The full number of parameters is number of points - the number of points corresponding to the reference

	
	//int diff_df = 2*(number_of_lineages-1); I think this is wrong.
	int diff_df = df - 2*(number_of_lineages-1);
	// Obtain second derivatives
	double dispersion_parameter = chi_square/diff_df;

	tsl::robin_map<string,double> gradient_s2_hash;  // Replicate independent
	// Obtain an approximation of the inverse of the Hessian.
	// We will just use the derivatives that are relevant to the lineage and assume no covariances between parameters across lineages.
	// To do so, we need the second derivatives with respect to s, all the f0, and the covariances between these parameters.
	// Recall that: d2l/ds2 = -sum(t^2 * N * f * (1-f))
	// Recall that: d2l/df02 = -sum(N * f * (1-f))
	// Recall that: d2l/dfs = -sum(t * N * f * (1-f))
	// All other covariances set to zero.
	// This is a symmetric arrowhead matrix. There exists a fast algorithm to compute this.
	// The modified Sherman-Morrison Inverse of arrowhead matrices can be used.
	// Let a matrix M = I + L + U
	// Where L are the strictly bottom elements of the matrix, and U are the strictly top bottom of the matrix.
	// The inverse of M is simply (I-L) * (I- 1/(1+alpha) * (U * (I-L) )
	// Alpha is -1 * U * L
	// In our case, we only care about the inverse of the s term, so the operation simplifies to 1/(1+alpha).
	// So we can easily find the inverse by dividing all the U elements by the d2l/ds2 term
	// Dividing all the L elements by the d2l/df02 terms.
	// Doing a dot product, and getting 1/(1+alpha).

	double tmp = 0;
	for(int i = 0;i < number_of_lineages;++i){
		vector<double> L;
		vector<double> U;
		double gradient_s = 0;
		double alpha = 0;
		// Is only computed from the data, not by replicates.
		int appearances = 0;
		for(int rep = 0;rep < number_of_replicates;++rep){
			double gradient_f = 0;
			double gradient_fs = 0;
			if(counts[rep].find(all_lineage_IDs[i]) != counts[rep].end()){
				int timepoints_size = timepoints[rep].size();
				++appearances;
				for(int t = 0;t < timepoints_size;++t){
					tmp = (double)total_reads[rep][t] *  expected[rep][all_lineage_IDs[i]][t] * (1-expected[rep][all_lineage_IDs[i]][t]);
					gradient_s -= timepoints[rep][t] * timepoints[rep][t] * tmp;										
					gradient_f -= tmp;
					gradient_fs -= timepoints[rep][t] * tmp;
				}
				// gradient_f are the diagonals
				// gradient_fs are the L and U.
				L.push_back(gradient_fs/gradient_f);
				U.push_back(gradient_fs);					
			}
			
		}

		for(int j = 0;j < appearances;++j){
			U[j] = U[j]/gradient_s;
			alpha += -1 * U[j] * L[j];
		}

		// inverse is 1 / (1+alpha) / gradient_s;
		gradient_s2_hash[all_lineage_IDs[i]] = 1 / (1+alpha) / gradient_s;

	}

	// Renormalize f0
	for(int rep = 0;rep < number_of_replicates;++rep){
		int timepoints_size = timepoints[rep].size();
		double sum_intercepts = 0.0;
		for(int i = 0;i < number_of_lineages;++i){
			if(counts[rep].find(all_lineage_IDs[i]) != counts[rep].end()){
				initial_frequencies[rep][all_lineage_IDs[i]] = exp(initial_frequencies[rep][all_lineage_IDs[i]]);
				sum_intercepts += initial_frequencies[rep][all_lineage_IDs[i]];
			}
		}
		for(int i = 0;i < number_of_lineages;++i){
			if(counts[rep].find(all_lineage_IDs[i]) != counts[rep].end()){
				initial_frequencies[rep][all_lineage_IDs[i]] /= sum_intercepts;		
			}
		}	
	}
	
	// Now have the standard error. Output the data!
	cout << "Lineage" << "	" << "s" << "	"<< "stderr(s)" << "	" << "f0_array" << endl;
	for(int i = 0;i < number_of_lineages;++i){
		cout << all_lineage_IDs[i] << "	" << selection_coefficients[all_lineage_IDs[i]] << "	";
		//if(i == number_of_lineages-1){
		//	cout << "0";
		//}
		//else{
			cout << sqrt(dispersion_parameter)*sqrt(-1*gradient_s2_hash[all_lineage_IDs[i]]);
		//}
		
		// Output the estimated f0.
		for(int rep = 0;rep < number_of_replicates;++rep){
			cout << "	";
			cout << initial_frequencies[rep][all_lineage_IDs[i]];
		}
		cout << endl;
	}
	cout << endl;
	time_1 = clock();
	walltime_1 = chrono::high_resolution_clock::now();
	chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(walltime_1 - walltime_0);

	// Output computational statistics
	cout << "Completed in " << iteration_count << " iterations. " << endl;
	cout << "Clock: " << diffclock(time_1,time_0)/1000 <<  "	" << "Wallclock: " << time_span.count() << endl;
	

	// Output some statistics.
	cout << "Data fit: " << endl;
	cout << "Log likelihood: " << log_likelihood << endl;
	cout << "Full likelihood: " << full_log_likelihood << endl;
	cout << "Scaled deviance: " << 2*(full_log_likelihood - log_likelihood)/dispersion_parameter << endl;
	cout << "df: " << diff_df << endl;
	cout << "Variance scaling factor " << dispersion_parameter << endl;
}