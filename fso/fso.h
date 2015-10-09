/*
 * fso.h
 *
 *  Created on: Dec 15, 2014
 */
#include <vector>
#ifndef FSO_H_
#define FSO_H_
//#define DEBUG
using namespace std;
class fso {
public:
	//fso(int dim, int n_individuals, int n_iterations, vector <double> parameters, int n_followed, vector< vector<double> > domain_, double max_vel_);
	void set_value(int dim, int n_individuals, int n_iterations, vector <double> parameters, int n_followed, vector< vector<double> > domain_, double max_vel_);
	virtual double compute_fitness_(vector<double>);
	// Init methods
	bool init();

	// Run methods
	bool compute_next_step();
	void debug();

private:

	bool init_matrix_interconnection();
	bool init_positions_();
	bool init_velocities_();
	bool init_global_best_values_();
	bool init_global_best_positions_();
	bool init_personal_best_values_();
	bool init_personal_best_positions_();
	bool compute_global_best_();
	bool compute_personal_best_();
	bool compute_velocities_();
	int control_domain();
	double compute_mean_velocity(int j, int i);
	vector< vector<unsigned int> > generate_mask_matrix(unsigned int n_one);
	unsigned int dim_;
	unsigned  n_individuals_;
	unsigned  n_iterations_;
	unsigned  n_followed_;
	double max_vel_;
	vector<double >fitness_;
	vector< vector<double> > domain_;
	vector<double> parameters_;
	vector< vector<double> > velocities_;
	vector< vector<double> > positions_;
	double global_best_value_;
	vector< double> global_best_position_;
	vector<double> personal_best_value_;
	vector< vector<double> > personal_best_position_;
	vector< vector<unsigned int> > matrix_interconnection_;

};

#endif /* FSO_H_ */
