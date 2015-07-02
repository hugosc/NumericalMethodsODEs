#include <limits>
#include <iostream>
#include <ctime>
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

typedef std::vector<double> state_type;

double max(double a, double b) {
	if (a>b) return a;
	else return b;
}

template <class Method, class System>
inline void benchmark_simulation (System system, state_type state, interval i, double step, std::pair<double,double> result) {
	Method method;
	clock_t init_time = clock();
	state_type solution = solve(method, system, state, i, step).back().second;
	clock_t end_time = clock();
	double duration = (static_cast<double>(end_time - init_time))/CLOCKS_PER_SEC;
	std::cout << "Simulation done in: " << duration << "\nfinal result: (" << solution[0] << "," << solution[1] << ")\nError: " << max(fabs(result.first - solution[0]), fabs(result.second - solution[1])) << "\n";
}

int main() {

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
	double step = 0.01;
	auto state2 = state;
	eval(runge_kutta_4th<state_type>(), system, state2, i, 0.00001);
	auto result = std::make_pair(state2[0], state2[1]);
	std::cout << "forward euler...\n----------------\n";
	benchmark_simulation<forward_euler<state_type>>(system, state, i, step, result);
	std::cout << "--------------\nbackward euler...\n---------------\n";
	benchmark_simulation<backward_euler<state_type>>(system, state, i, step, result);
	std::cout << "--------------\nmodified euler...\n---------------\n";
	benchmark_simulation<modified_euler<state_type>>(system, state, i, step, result);
	std::cout << "--------------\nrunge-kutta 2nd order...\n---------------\n";
	benchmark_simulation<runge_kutta_2nd<state_type>>(system, state, i, step, result);
	std::cout << "--------------\nrunge-kutta 3rd order...\n---------------\n";
	benchmark_simulation<runge_kutta_3rd<state_type>>(system, state, i, step, result);
	std::cout << "--------------\nrunge-kutta 4th order...\n---------------\n";
	benchmark_simulation<runge_kutta_4th<state_type>>(system, state, i, step, result);
	std::cout << "--------------\npredictor-corrector adams 3rd order...\n---------------\n";
	benchmark_simulation<predictor_corrector_3rd<state_type>>(system, state, i, step, result);
	std::cout << "--------------\npredictor-corrector adams 4th order...\n---------------\n";
	benchmark_simulation<predictor_corrector_4th<state_type>>(system, state, i, step, result);
	return 0;
}
