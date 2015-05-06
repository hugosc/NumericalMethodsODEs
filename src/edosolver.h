#include <functional>

/**
* @file edosolver.h
* @brief file containing the interface for the solving methods,
* such as Euler methods or Runge-Kutta methods
*/

/**
* @brief A lambda function that receives and returns
* a double is used to represent a real function
*/
typedef std::function<double (double)> realfunc;


/**
 * @brief interface for the intial value
 * problem solvers
 * @todo properly define the class
 */
class edo_solver {
	public:
		virtual  edo_solver() = 0;
		virtual ~edo_solver() = 0;
	protected:	
}