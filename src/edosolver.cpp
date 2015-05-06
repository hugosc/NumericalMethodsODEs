#include "edosolver.h"

/**
 * @file edosolver.cpp
 * @brief implementation of the edo_solver class
 */

void edo_solver::set_step_size(double sz) {
	_step_size = sz;
}

double edo_solver::step_size() {
	return _step_size;
}


std::vector<std::pair<double,double>> edo_solver::solve (interval time_interval) {
	double current_time = time_interval.first;
	int number_of_intervals = 1 + static_cast<int>((time_interval.second-time_interval.first)/_step_size);
	std::vector<std::pair<double,double>> solutions;
	solutions.reserve(number_of_intervals);
	while(time_interval.second - current_time > epsilon) {
		solutions.push_back(std::make_pair(current_time,approx_point(current_time)));
		current_time += _step_size;
	}
	return solutions;
}