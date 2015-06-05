#ifndef TYPES_H
#define TYPES_H

#include <cmath>

struct interval {
	double first, second;
	static constexpr double epsilon = std::numeric_limits<double>::epsilon();
	static double eq_eps(const double a, const double b) {
		return fabs(a-b) < epsilon;
	}
	interval(double f, double s) : first(f), second(s) {}
	bool contains (const double a) {
		return (a>=first && a<=second);
	}
};

#endif //TYPES_H
