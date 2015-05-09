#include "odesolver.h"

/**
 * @file odesolver.cpp
 * @brief implementation of the ode_solver class
 */

void ode_solver::set_step_size(double sz) {
	_step_size = sz;
}

double ode_solver::step_size() {
	return _step_size;
}

double ode_solver::eval (double x) {
	double current_val  = _problem.initial_vals.front();
	for (double i=_problem.initial_time ; !equals_by_eps(i,x) && i<x ; i+=_step_size) {
		current_val = approx_point(i,current_val);
	}
	return current_val;
}

std::vector<std::pair<double,double>> ode_solver::solve (interval time_interval) {
	double current_time = time_interval.first;
	int number_of_intervals = 1 + static_cast<int>((time_interval.second-time_interval.first)/_step_size);
	std::vector<std::pair<double,double>> solutions;
	solutions.reserve(number_of_intervals);
	double current_val = eval(current_time);
	while(!equals_by_eps(current_time,time_interval.second) && current_time<time_interval.second) {
		solutions.push_back(std::make_pair(current_time,current_val));
		current_time += _step_size;
		current_val = approx_point(current_time,current_val);
	}
	return solutions;
}

double forward_euler::approx_point(double t, double last_point) {
	return last_point + _step_size*_problem.equation(last_point,t);
}