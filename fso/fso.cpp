/*
 * fso.cpp
 *
 *  Created on: Dec 15, 2014
 *      Author: pulcini
 */
#include <vector>
#include <iostream>
#include "fso.h"
#include <algorithm>
#include "assert.h"
#include <iostream>
#include <fstream>
#include <climits>
#include <math.h>
#include <omp.h>
void fso::set_value(int dim, int n_individuals, int n_iterations, vector<double> parameters, int n_followed,vector< vector<double> > domain, double max_vel) {
	dim_ = dim;
	n_individuals_ = n_individuals;
	n_iterations_ = n_iterations;
	parameters_ = parameters;
	n_followed_ = n_followed;
	domain_ = domain;
	max_vel_ = max_vel;
	// init
	vector<double> init_;
	////cout << "DEBUG 1" << endl;
	for(unsigned int i = 0; i < dim_; i++) {
	    /* std:://cout << a[i]; ... */
		init_.push_back(0.0);
	}
	//cout << "DEBUG 2" << endl;
	for(unsigned int i = 0; i < n_individuals_; i++) {
	    /* std:://cout << a[i]; ... */
		velocities_.push_back(init_);
		positions_.push_back(init_);
	}
	//cout << "DEBUG 3" << endl;
	vector<int>  matrix_;
	for(unsigned int i = 0; i < n_followed_; i++) {
	    /* std:://cout << a[i]; ... */
		matrix_.push_back(0);
	}
	//cout << "DEBUG 4" << endl;
	vector<unsigned int> matrix_int;
	//cout << "DEBUG 5" << endl;
	global_best_value_ = 1000000000;
	//cout << "DEBUG 6" << endl;
}


bool fso::init() {
	// TODO Auto-generated destructor stub
	//cout << "DEBUG 7" << endl;
	if(!init_matrix_interconnection())
		return false;
	//cout << "DEBUG 8" << endl;
	if(!init_positions_())
		return false;
	//cout << "DEBUG 9" << endl;
	if(!init_velocities_())
		return false;
	//cout << "DEBUG 10" << endl;
	if(!init_global_best_values_())
		return false;
	//cout << "DEBUG 11" << endl;
	if(!init_global_best_positions_())
		return false;
	//cout << "DEBUG 12" << endl;
	if(!init_personal_best_values_())
		return false;
	//cout << "DEBUG 13" << endl;
	if(!init_personal_best_positions_())
		return false;
	//cout << "DEBUG 14" << endl;

	for(unsigned int i = 0; i < n_individuals_; ++i){
		fitness_.push_back(compute_fitness_(positions_.at(i)));
	}

	return true;

}
bool fso::init_matrix_interconnection(){
	for( unsigned int i = 0; i < n_individuals_;++i){
		vector<unsigned int> mask;
		for(unsigned f = 0; f < n_individuals_; ++f){
			if(f!=i){ mask.push_back(f); } };
		random_shuffle(mask.begin(),mask.end());
		matrix_interconnection_.push_back(mask);
	}

	return true;
}
bool fso::init_positions_(){

	for (unsigned int j=0; j < n_individuals_; ++j){
		for(unsigned int i=0; i < dim_; ++i){
			//cout << ((double) rand() / (RAND_MAX)) * domain_.at(1).at(i) << endl;
			positions_.at(j)[i] = ((double) rand() / (RAND_MAX)) * domain_.at(1).at(i);
		}
		for(unsigned int i = 0; i < n_individuals_; ++i){
#ifdef DEBUG
				//cout << "VEL" << endl;
				//cout << positions_.at(i)[0] << endl;
#endif
		}
	}
	return true;
}
bool fso::init_velocities_(){

	for (unsigned int j=0; j < n_individuals_; ++j){
		for (unsigned int i=0; i < dim_; ++i){
			velocities_.at(j)[i] = (((double) rand() / (RAND_MAX)));
		}
	}

	for(unsigned int i = 0; i < n_individuals_; ++i){
		//cout << "VEL" << endl;
		//cout << velocities_.at(i)[0] << endl;
	}
	return true;
}
bool fso::init_global_best_values_(){
	global_best_value_ = 0.0;
	return true;
}
bool fso::init_global_best_positions_(){
	vector<double> internal (dim_, 10.0);
	global_best_position_ = internal;
	return true;
}

bool fso::init_personal_best_values_(){
	vector<double> internal (n_individuals_, LONG_MAX);
	personal_best_value_ = internal;
	return true;
}
bool fso::init_personal_best_positions_(){
	vector<double> internal (dim_, 0.0);
	for (unsigned int j = 0; j < n_individuals_; ++j){
		personal_best_position_.push_back(internal);
	}
	return true;
}

bool fso::compute_global_best_(){
	int min_index = min_element(personal_best_value_.begin(), personal_best_value_.end())
			- personal_best_value_.begin();
	global_best_value_ = personal_best_value_.at(min_index);
	global_best_position_ = personal_best_position_.at(min_index);
	return true;
}

bool fso::compute_personal_best_(){
	//vector<double> fitness;

	unsigned int i;
	for(i= 0; i < n_individuals_; ++i){
		if(fitness_.at(i) < personal_best_value_.at(i)){
			personal_best_value_.at(i) = fitness_.at(i);
			personal_best_position_.at(i) = positions_.at(i);
		}
	}
	return true;
}

bool fso::compute_velocities_(){
	vector< vector<double> > new_velocities;
	vector< vector<double> > new_position = positions_;
	unsigned int i,j;
//        #pragma omp parallel shared(new_velocities,new_position) private(j,i)
  //      #pragma omp for
	for(j = 0; j < n_individuals_ ; ++j){
		vector<double> par;
		for(i = 0; i < dim_; ++i){
			////cout << "J= " << j << endl;
			////cout << "I= " << i << endl;
			double par;
			////cout << "PAR X = " << parameters_.at(0) * velocities_.at(j)[i] << endl;
			////cout << "PAR XX = " << velocities_.at(j)[i] << endl;
			par = - parameters_.at(0) * velocities_.at(j)[i]+
					+ parameters_.at(1) * (global_best_position_.at(i) - positions_.at(j)[i]) +
					+ parameters_.at(2) * (personal_best_position_.at(j)[i] - positions_.at(j)[i]) +
					+ compute_mean_velocity(j,i);
			////cout << "PAR= " << par << endl;
			////cout << "DEBUG 33" << endl;
			if( par > max_vel_){
				if(((double) rand() / (RAND_MAX)) > 0.5){
					par = ((double) rand() / (RAND_MAX));
			}else{
				par = -1 * ((double) rand() / (RAND_MAX));
			}
			}
			new_position.at(j).at(i) = new_position.at(j).at(i) + par;
			if(new_position.at(j).at(i) > domain_.at(1).at(i) || new_position.at(j).at(i) < domain_.at(0).at(i)){
				////cout << par << endl;
				par = -par;
				////cout << par << endl;
				////cout << positions_.at(j).at(i) << endl;
				////cout << new_position.at(j).at(i) << endl;
				positions_.at(j).at(i) = positions_.at(j).at(i) + par;
				////cout << "BìNEGATIVO" << endl;
			}else{
				positions_.at(j).at(i) = positions_.at(j).at(i) + par;
			}


			velocities_.at(j).at(i) = par;
	}
		////cout << "DEBUG 33 fuori" << endl;

}
	return true;
}


double fso::compute_mean_velocity(int j, int i){
	double ratio = (double)1/(double)n_followed_;
	double sum = 0;
	for(unsigned int k = 0; k < n_followed_; ++k ){
		////cout << velocities_.at(matrix_interconnection_.at(j).at(k)).at(i) << endl;
		sum = sum + velocities_.at(matrix_interconnection_.at(j).at(k)).at(i);
	}
	sum = sum * ratio;
	////cout << sum << endl;
	////cout << "fuori" << endl;
	return sum;
}

bool fso::compute_next_step(){
	compute_velocities_();
	////cout << "DEBUG 31" << endl;

	for(unsigned int i = 0; i < n_individuals_; ++i){
		fitness_.at(i) = compute_fitness_(positions_.at(i));
		/*//cout << compute_fitness_(positions_.at(i)) << endl;
		//cout << fitness_.at(i) << endl;*/
	}
	compute_personal_best_();
	compute_global_best_();
	return true;
}

double fso::compute_fitness_(vector<double> vec){
	//Ackley
	double f;
	//f =-20*exp(-0.2* sqrt(0.5*( pow(vec.at(0),2) + pow(vec.at(0),2) )) ) -
	//		exp(0.5*( cos(2*3.14*vec.at(0)+cos(2*3.14*vec.at(1)))));
	/*
	 * Beale's function
	 */
	double par1,par2;
	par1 = pow(vec.at(0)*vec.at(1),2);
	par2 = pow(vec.at(0)*vec.at(1),3);
	f = pow((1.5 - vec.at(0) + vec.at(0)*vec.at(1)),2) + pow((2.25-vec.at(0)+par1),2)+
			pow((2.25-vec.at(0)+par2),2);
	/*
	 *
	 */
	return f;
	//return sqrt(pow(vec.at(0),2)+pow(vec.at(1),2)+3);
	//return pow(vec.at(0),3)*log(vec.at(1))+3.4;
}

void fso::debug(){
	  ofstream position1,velocity, global_best,gbp;

	  position1.open ("position1.txt",ios::app);

	  gbp.open ("gbp.txt",ios::app);
	  velocity.open ("velocity.txt",ios::app);
	  global_best.open("global_best.txt",ios::app);

	  for(unsigned int i=0; i < n_individuals_; ++i){
		  position1 << positions_.at(i)[0] << " " << positions_.at(i)[1] << " " ;
	  }
	  position1 << endl;

	  gbp <<  global_best_position_.at(0) << " " <<  global_best_position_.at(1) << " ";
	  gbp << global_best_value_ << endl;

	  for(unsigned int i=0; i < n_individuals_; ++i){
		  velocity << velocities_.at(i)[0] << " " << velocities_.at(i)[1] << " " ;
	  }
		velocity << endl;
		global_best << global_best_value_ << endl;

		global_best.close();
		velocity.close();
		gbp.close();
}
