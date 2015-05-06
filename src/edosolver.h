#ifndef EDOSOLVER_H
#define EDOSOLVER_H

#include <vector>
#include <utility>
#include <limits>
#include "types.h"

/**
* @file edosolver.h
* @brief file containing the interface for the solving methods,
* such as Euler methods or Runge-Kutta methods
*/

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
 * @brief interface for the intial value
 * problem solvers
 * @todo Decide if more member functions are to be added
 */
class edo_solver {
	public:
		constexpr static const double epsilon = std::numeric_limits<double>::epsilon(); /**<constant for machine epsilon used for comparing doubles*/
		/**
		 * @brief Class constructor
		 * @param problem problem to be solved
		 * @param sz step size for the solution
		 */
		edo_solver(initval_problem problem, double sz) : _problem(problem), _step_size(sz) {}
		/**
		 * @brief function to solve the init val problem
		 * @param time_interval time interval for the solution
		 * @return vector of pairs (t,y(t)) 
		 * y being the solution to the EDO
		 * time starts at beggining of interval and increments by step size
		 * up to the end of the interval, and at each iteration, it makes a call to
		 * approx_point and pushes its return value to the vector that will be returned
		 */
		std::vector<std::pair<double,double>> solve (interval time_interval);
		
		/** * @brief setter for _step_size */
		virtual void set_step_size(double sz);
		/** 
		* @brief getter for _step_size
		* @return _step_size
		*/
		double step_size();
	protected:
		double _step_size; /**<how much the time is incremented in each interation */
		std::vector<std::pair<double,double>> _solutions; /**<stores the solutions to the problem */
		initval_problem _problem; /**<problem to be solved */

		/**
		 * @brief finds an approximation for the solution function
		 * applied to t
		 * @param t specific time
		 * @return value of solution function applied to t
		 */
		virtual double approx_point(double t) = 0;
};

#endif //EDOSOLVER_H
