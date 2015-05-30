#ifndef TYPES_H_
#define TYPES_H_

#include <functional>
#include <vector>
#include <utility>
#include <initializer_list>

/**@file types.h 
* @brief definition of important types, such as 
* real functions(realfunc) and Initial Value Problem(initval_problem)
*/

/**
* @brief A lambda function that receives and returns
* a double is used to represent a real function
*/
typedef std::function<double (double)> realfunc;

/**
 * @brief explicit differential equation
 */
typedef std::function<double (double,double)> diff_equation;

/**
 * @brief struct to represent a closed interval
 */
struct interval {
	double first; /**<lim inf of interval */
	double second; /**<lim sup of interval */
	/**
	 * @brief creates interval, default interval is [0,0].
	 * @param f sets member first to its value
	 * @param s sets member second to its value
	 */
	interval (double f=0.0, double s=0.0) : first(f), second(s) {}
};

/**
 * @brief class to represent an initial value problem
 * @todo Decide if member order should be removed and allow only first order differential
 * equations
 */
class initval_problem {
	public:
		const diff_equation equation; /**<Ordinary differential equation for the problem*/
		const std::vector<double> initial_vals; /**<initial values for the problem*/
		const double initial_time; /**<initial time t0 for the problem*/
		const int order; /**<Order of the differential equation*/
		/** @brief constructor */
		initval_problem(diff_equation e,std::vector<double> iv, double it) : equation(e), initial_vals(iv), initial_time(it), order(iv.size()) {}
		/** @brief constructor with initializer list */
		initval_problem(diff_equation e,std::initializer_list<double> iv, double it) : equation(e), initial_vals(iv), initial_time(it), order(iv.size()) {}
};

#endif //TYPES_H_