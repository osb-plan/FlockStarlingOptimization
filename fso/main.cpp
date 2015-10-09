#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <time.h>
#include <string>
#include <stdlib.h>
//#include "myfso.h"
#include "fso.h"
#include <omp.h>

using namespace std;
int main(int argc, char *argv[]){



	vector<int> dentro, dopo;
	vector<vector<int> > mio;
	dentro.push_back(1);
	dentro.push_back(2);
	dentro.push_back(3);
	dopo.push_back(1);
	dopo.push_back(2);
	dopo.push_back(3);

	mio.push_back(dentro);
	mio.push_back(dopo);
	//cout << mio.at(0)[0] << endl;
	//cout << mio[0][2] << endl;
	//cout << mio.at(0)[0] << endl;
	//cout << mio.at(0)[0] << endl;
	//return 0;
	/* Parameters
	 *
	 */
	int dim = 2;
	std::string::size_type sz;

	int n_individuals = stoi(argv[1],&sz);
	int n_iterations = 1000;
	int n_followed = 7;
	vector<double > parameters;

	parameters.reserve(dim);
	parameters.push_back(0.25);
	parameters.push_back(0.005);
	parameters.push_back(0.0005);

	vector< vector<double> > domain;
	//domain.reserve(dim);
	vector<double> min_dom (dim, -10.0);
	vector<double> max_dom (dim, 10.0);
	for(int j = 0; j < dim; ++j){
	domain.push_back(min_dom);
	domain.push_back(max_dom);
		/*domain.at(0).at(j) = -10;
		domain.at(1).at(j) = 10;*/
	}
	//cout << "pippo" << endl;
	double max_vel = 0.5;
	fso FS;
	FS.set_value(dim,n_individuals,n_iterations,parameters,n_followed,domain, max_vel);
	FS.init();
    //  Start Timers
    clock_t wall0,wall1;
    clock_t cpu0,cpu1;
    //  Start Timers
    cpu0 = clock();

    //omp_set_num_threads(4);
    int k;
//    #pragma omp parallel shared(FS) private(k)
//    #pragma omp for
    for(k = 0; k < n_iterations; ++k){


		FS.compute_next_step();


		//cout << k << endl;
#ifdef DEBUG
		FS.debug();
#endif
	}
    cpu1  = clock();
    cout <<  ((float)(cpu1 -cpu0))/CLOCKS_PER_SEC<< endl;

	//myfso FSM (dim,n_individuals,n_iterations,parameters,n_followed,domain, max_vel);
	/*myfso FSM;
	FSM.set_value(dim,n_individuals,n_iterations,parameters,n_followed,domain, max_vel);
	vector<double>  a;
	cout << FSM.compute_fitness_(a) << endl;
	cout << FS.compute_fitness_(a) << endl;*/

}
