#include <limits>
#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <array>
#include "types.h"
#include "num_methods.h"

//diretiva de compilacao : g++ exemplo.cpp -o exemplo -std=c++11 -O3

static const double k1 = 100;
static const double k2 = 100;
static const double m1 = 1;
static const double m2 = 1;
static const double g = 9.81;
static const double b1 = 10;
static const double b2 = 10;

int main() {
	typedef std::vector<double> state_type;

	auto system = [=](const state_type& x, double t) {
		state_type dxdt(4);
		dxdt[0] = x[2];
		dxdt[1] = x[3];
		dxdt[2] = g - ((b1 + b2)/m1)*x[2] - ((k1 + k2)/m1)*x[0] + (b2/m1)*x[3] + (k2/m1)*x[1]; 
		dxdt[3] = g + (b2/m2)*x[2] + (k2/m2)*x[0] - (b2/m2)*x[3] - (k2/m2)*x[1];
		return dxdt;
	};
	interval i(0,20);
	state_type state{{0,1,0,0}};
	predictor_corrector_3rd<state_type> method;
//	runge_kutta2nd<state_type> method;
	auto simulation = solve(method, lambda, state, i, 0.001);
	for (auto pair : v) {
		std::cout << "(" << pair.second[0] << "," << pair.second[1] << ")\n";
	}
	return 0;
}
