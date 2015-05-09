#include "odesolver.h"
#include <iostream>

int main() {
	initval_problem problem([](double y, double x) {
		return sqrt(y) + x;
	},{1.0},0);
	auto method = forward_euler(problem,0.000001);
	std::cout << method.solve(interval(0,20)).back().second << std::endl;
	return 0;
}