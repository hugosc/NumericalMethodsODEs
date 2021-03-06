#ifndef ODESOLVER_H
#define ODESOLVER_H

#include <vector>
#include <utility>
#include <limits>
#include <cmath>
#include "types.h"

/**
* @file odesolver.h
* @brief file containing the interface for the solving methods,
* such as Euler methods or Runge-Kutta methods
*/

/**
 * @brief interface for the intial value
 * problem solvers
 * @todo Decide if more member functions are to be added
 */
class ode_solver {
	public:
		constexpr static const double epsilon = std::numeric_limits<double>::epsilon(); /**<constant for machine epsilon used for comparing doubles*/
		/**
		 * @brief Class constructor
		 * @param problem problem to be solved
		 * @param sz step size for the solution
		 */
		ode_solver(initval_problem problem, double sz) : _problem(problem), _step_size(sz) {}
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

		/**
		 * @brief returns a lambda that is the approximate solution to the problem
		 * @return real function that is the solution to the problem
		 */
		realfunc solution_func();
		
		/** * @brief setter for _step_size */
		virtual void set_step_size(double sz);
		/** 
		* @brief getter for _step_size
		* @return _step_size
		*/
		double step_size();

		/**
		 * @brief operator overload for eval(x)
		 * @param x value to apply in the solution function
		 * @return eval(x)
		 */
		double operator() (double x) { return eval(x); }

		/** @brief checks if two double are within the computer epsilon of difference */
		static bool equals_by_eps(const double x, const double y) {
			return fabs(x-y) < epsilon;
		}
	protected:
		double _step_size; /**<how much the time is incremented in each interation */
		initval_problem _problem; /**<problem to be solved */
		/**
		 * @brief gives an approximation to y(x), y being the solution to the problem
		 * @param x value to apply in the solution function
		 * @return approximation for y(x)
		 */
		double eval(double x);	
		/**
		 * @brief finds an approximation for the solution function
		 * applied to t
		 * @param t specific time
		 * @param last_point last obtained appximation
		 * @return value of solution function applied to t
		 */
		virtual double approx_point(double t, double last_point) = 0;
};

/**
 * @brief class that inherits the ode_solver interface
 * and implements the forward euler method
 */
class forward_euler : public ode_solver {
	public:
		/** @brief constructor calls its parent constructor */
		forward_euler(initval_problem problem, double sz) : ode_solver(problem,sz) {}
		/**
		 * @brief gives an approximation to y(x), y being the solution to the problem
		 * @param x value to apply in the solution function
		 * @return approximation for y(x)
		 */
		double operator() (double x) { eval(x); } 
	protected:
		/**
		 * @brief finds an approximation for the solution function
		 * applied to t
		 * @param t specific time
		 * @param last_point last obtained approximation
		 * @return value of solution function applied to t
		 */
		double approx_point(double t,double last_point);
};

/**
*@brief class that inherits the ode_solver interface and implements the second order Runge-Kutta method
* 
*/
class secondorder_rungekutta : public ode_solver
{
	public:
		/** @brief constructor calls ode_solver()*/
		secondorder_rungekutta(initval_problem problem, double sz) : ode_solver(problem,sz){}
		/** @brief function object that evaluate the approximation of y(x)
		 *  @param x value to apply in the solution function
		 *  @return  approximation for y(x)*/
		double operator() (double x) {eval(x);}
	protected:
		
		/** @brief finds an approximation for the solution function applied to t
		 * @param t specific time
		 * @param last_point last obtained approximation
		 */
		double approx_point(double t, double last_point);
};

#endif //ODESOLVER_H
