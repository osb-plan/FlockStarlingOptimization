#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <time.h>
#include <string>
//#include "myfso.h"
#include "fso.h"
#include <omp.h>

using namespace std;
int main(int argc, char *argv[]){

	/* Parameters
	 *
	 */
	int dim = 2;
	std::string::size_type sz;
	// CTRL input data
	if(argc<=1){
                fprintf(stderr,"Error you have to provide number of individuals!\n");
                exit(EXIT_FAILURE);
        }

	int n_individuals = stoi(argv[1],&sz);
	int n_iterations = 5000;
	int n_followed = 7;


	if( n_individuals <= n_followed){
		fprintf(stderr,"Error nr. of individuals must be greater than followed birds\n");
		exit(EXIT_FAILURE);
	}

	vector<double > parameters;

	parameters.reserve(dim);
	parameters.push_back(1);
	parameters.push_back(0.1);
	parameters.push_back(0.05);

	vector< vector<double> > domain;
	//domain.reserve(dim);
	vector<double> min_dom (dim, -40.0);
	vector<double> max_dom (dim, 40.0);
	for(int j = 0; j < dim; ++j){
	domain.push_back(min_dom);
	domain.push_back(max_dom);
		/*domain.at(0).at(j) = -10;
		domain.at(1).at(j) = 10;*/
	}
	//cout << "pippo" << endl;
	double max_vel = 1;
	fso FS;
	FS.set_value(dim,n_individuals,n_iterations,parameters,n_followed,domain, max_vel);
	FS.init();
    	int k;
    	for(k = 0; k < n_iterations; ++k){


		FS.compute_next_step();


	//cout << k << endl;
	#ifdef DEBUG
			FS.debug();
	#endif
	}
	
	//myfso FSM (dim,n_individuals,n_iterations,parameters,n_followed,domain, max_vel);
	/*myfso FSM;
	FSM.set_value(dim,n_individuals,n_iterations,parameters,n_followed,domain, max_vel);
	vector<double>  a;
	cout << FSM.compute_fitness_(a) << endl;
	cout << FS.compute_fitness_(a) << endl;*/

}
