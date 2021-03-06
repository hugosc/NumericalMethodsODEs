#ifndef NUM_METHODS_H
#define NUM_METHODS_H

#include "types.h"
#include <vector>
#include <limits>
#include <iostream>

template <class State>
void print_state(const State& state) {
	std::cout << "[ ";
	for (auto v : state)
		std::cout << v << " ";
	std::cout << "]\n";
}

template <class State>
double max_diff(const State& a, const State& b) {
	double d = -1;
	for (int i = 0; i < a.size(); ++i) {
		double diff = fabs(a[i] - b[i]);
		if (diff > d) d = diff;
	}
	return d;
}

template <class Method>
struct method_traits {
	static const bool is_method = false;
};

template <class State>
class forward_euler {
	public:
		template <class System>
		void approx_point(System system, State& state, double curr_time, double step) {
			auto dxdt = system(state, curr_time);
			for (int i=0 ; i<state.size(); i++)
				state[i] += step * dxdt[i];
		}

		inline void cleanup() {}
};

template <class State>
struct method_traits<forward_euler<State>> {
	static const bool is_method = true;
};

template <class State>
class backward_euler {
	public:
		template <class System>
		void approx_point(System system, State& state, double curr_time, double step) {
			auto curr_state = state;
			State prev_state;
			do {
				prev_state = curr_state;
				auto dxdt = system(curr_state, curr_time + step);
				for (int i = 0; i < state.size(); ++i)
					curr_state[i] = state[i] + step*dxdt[i];
			} while (max_diff(curr_state, prev_state) > step);
			state = curr_state;
		}

		inline void cleanup() {}
};

template <class State>
struct method_traits<backward_euler<State>> {
	static const bool is_method = true;
};


template <class State>
class modified_euler {
	public:
		template <class System>
		void approx_point(System system, State& state, double curr_time, double step) {
			auto curr_state = state;
			State prev_state;
			do {
				prev_state = curr_state;
				auto dxdt1 = system(curr_state, curr_time + step);
				auto dxdt2 = system(state, curr_time);
				for (int i = 0; i < state.size(); ++i)
					curr_state[i] = state[i] + (step/2.0)*(dxdt1[i] + dxdt2[i]);
			} while (max_diff(curr_state, prev_state) > step);
			state = curr_state;
		}

		inline void cleanup() {}
};

template <class State>
struct method_traits<modified_euler<State>> {
	static const bool is_method = true;
};

template <class State>
class predictor_corrector_3rd {
	public:
		template <class System>
		void approx_point(System system, State& state, double curr_time, double step) {
			switch(iter) {
				case 0:
					y0_ = system(state, curr_time);
					iter++;
					break;
				case 1:
					y1_ = system(state, curr_time);
					for (int i=0; i<state.size(); i++) {
						state[i] += step * y1_[i]; 
					}
					iter++;
					break;
				case 2:
					y2_ = system(state, curr_time);
					for (int i=0; i<state.size(); i++) {
						state[i] += step * y2_[i];
					}
					iter++;
					break;
				case 3:
					auto y2 = state;
					for (int i=0; i<state.size(); i++) {
						state[i] = y2[i] + (step/12)*(23*y2_[i] - 16*y1_[i] + 5*y0_[i]);
					}
					y0_ = y1_;
					y1_ = y2_;
					y2_ = system(state,curr_time);
					for (int i=0; i<state.size(); i++) {
						state[i] = y2[i] + (step/12)*(5*y2_[i] + 8*y1_[i] -y0_[i]);
					}
					break;
			}
		}
		inline void cleanup(){iter=0;}
	private:
		State y0_, y1_, y2_;
		int iter = 0;
};

template <class State>
struct method_traits<predictor_corrector_3rd<State>> {
	static const bool is_method = true;
};

template <class State>
class predictor_corrector_4th {
	public:
		template <class System>
		void approx_point(System system, State& state, double curr_time, double step) {
			switch(iter) {
				case 0:
					y0_ = system(state, curr_time);
					iter++;
					break;
				case 1:
					y1_ = system(state, curr_time);
					for (int i=0; i<state.size(); i++) {
						state[i] += step * y1_[i];
					}
					iter++;
					break;
				case 2:
					y2_ = system(state, curr_time);
					for (int i=0; i<state.size(); i++) {
						state[i] += step * y2_[i];
					}
					iter++;
					break;
                case 3:
                    y3_ = system(state, curr_time);
					for (int i=0; i<state.size(); i++) {
						state[i] += step * y3_[i];
					}
					iter++;
					break;

				case 4:
					auto y3 = state;
					for (int i=0; i<state.size(); i++) {
						state[i] = y3[i] + (step/24)*(55*y3_[i] - 59*y2_[i] + 37*y1_[i] - 9*y0_[i]);
					}
					y0_ = y1_;
					y1_ = y2_;
					y2_ = y3_;
					y3_ = system(state,curr_time);
					for (int i=0; i<state.size(); i++) {
						state[i] = y3[i] + (step/24)*(9*y3_[i] + 19*y2_[i] - 5*y1_[i] + y0_[i]);
					}
					break;
			}
		}
		inline void cleanup(){iter=0;}
	private:
		State y0_, y1_, y2_, y3_;
		int iter = 0;
};

template <class State>
struct method_traits<predictor_corrector_4th<State>> {
	static const bool is_method = true;
};

template <class State>
class runge_kutta_2nd {
	public:
		template <class System>
		void approx_point(System system, State& state, double curr_time, double step) {

			auto k1 = system(state, curr_time);
			for (int i = 0; i < state.size(); i++)
				k1[i] *= step;

			auto temp = k1;

			for (int i = 0; i < state.size(); i++)
				temp[i] += state[i];

			auto k2 = system(temp, curr_time + step);

			for (int i = 0; i < state.size(); i++)
				k2[i] *= step;

			for (int i = 0; i < state.size(); i++)
				state[i] += k1[i]/2.0 + k2[i]/2.0;
		}

		inline void cleanup() {}
};

template <class State>
struct method_traits<runge_kutta_2nd<State>> {
	static const bool is_method = true;
};

template <class State>
class runge_kutta_3rd {
	public:
		template <class System>
		void approx_point(System system, State& state, double curr_time, double step) {

			auto k1 = system(state, curr_time);
			for (int i = 0; i < state.size(); ++i)
				k1[i] *= step;

			auto temp = k1;
			for (int i = 0; i < state.size(); ++i)
				temp[i] = state[i] + temp[i]/2.0;

			auto k2 = system(temp, curr_time + step/2.0);
			for (int i = 0; i < state.size(); ++i)
				k2[i] *= step;

			for (int i = 0; i < state.size(); ++i)
				temp[i] = state[i] - k1[i] + 2*k2[i];

			auto k3 = system(temp, curr_time + step);
			for (int i = 0; i < state.size(); ++i)
				k3[i] *= step;

			for (int i = 0; i < state.size(); ++i)
				state[i] += (k1[i] + 4*k2[i] + k3[i])/6.0; 
		}
		inline void cleanup() {}
};

template <class State>
struct method_traits<runge_kutta_3rd<State>> {
	static const bool is_method = true;
};

template <class State>
class runge_kutta_4th {
	public:
		template <class System>
		void approx_point(System system, State& state, double curr_time, double step) {

			auto k1 = system(state, curr_time);
			for (int i = 0; i < state.size(); ++i)
				k1[i] *= step;

			auto temp = k1;
			for (int i = 0; i < state.size(); ++i)
				temp[i] = state[i] + temp[i]/3.0;

			auto k2 = system(temp, curr_time + step/3.0);
			for (int i = 0; i < state.size(); ++i)
				k2[i] *= step;

			for (int i = 0; i < state.size(); ++i)
				temp[i] += k2[i]/3.0;

			auto k3 = system(temp, curr_time + (2.0/3.0)*step);
			for (int i = 0; i < state.size(); ++i)
				k3[i] *= step;

			auto k4 = state;
			for (int i = 0; i < state.size(); ++i)
				k4[i] += (k1[i] - k2[i] + k3[i]);
			k4 = system(k4, curr_time + step);
			for (int i = 0; i < state.size(); ++i)
				k4[i] *= step;

			for (int i = 0; i < state.size(); ++i)
				state[i] += (k1[i] + 3*k2[i] + 3*k3[i] + k4[i])/8.0;
		}
		inline void cleanup() {}
};

template <class State>
struct method_traits<runge_kutta_4th<State>> {
	static const bool is_method = true;
};

template <class State, class System, class Method>
std::vector<std::pair<double,State>> solve(Method method, System system, State state, interval i, double step) {
	static_assert(method_traits<Method>::is_method, "Not a valid method");
	std::vector<std::pair<double,State>> approximations;
	approximations.reserve(int((i.second-i.first)/step)+1);
	auto curr_time = i.first;
	method.cleanup();
	approximations.push_back(std::make_pair(curr_time, state));
 	while (i.contains(curr_time)) {
		method.approx_point(system, state, curr_time, step);
		curr_time += step;
		approximations.push_back(std::make_pair(curr_time, state));
	}
	return approximations;
}

template <class State, class System, class Method>
void eval(Method method, System system, State& state, interval i, double step) {
	static_assert(method_traits<Method>::is_method, "Not a valid method");
	auto curr_time = i.first;
	method.cleanup();
	while (i.contains(curr_time)) {
		method.approx_point(system, state, curr_time, step);
		curr_time += step;
	}
}

#endif //NUM_METHODS_H
