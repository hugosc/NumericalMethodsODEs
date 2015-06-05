#ifndef TYPES_H
#define TYPES_H

#include <cmath>

/**@file types.h 
* @brief definition of important types, such as interval
*/

/**
 * @brief struct to represent a closed interval
 */
struct interval {
	double first, second;
	static constexpr double epsilon = std::numeric_limits<double>::epsilon();
	static double eq_eps(const double a, const double b) {
		return fabs(a-b) < epsilon;
	}
	 /**
	 * @brief creates interval, default interval is [0,0].
	 * @param f sets member first to its value
	 * @param s sets member second to its value
	 */
	interval(double f, double s) : first(f), second(s) {}
	bool contains (const double a) {
		return (a>=first && a<=second);
	}
};

#endif //TYPES_H
