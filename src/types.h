#ifndef TYPES_H_
#define TYPES_H_

#include <functional>

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
 * @brief class to represent an initial value problem
 * @todo properly define it
 */
class initval_problem {

};

#endif //TYPES_H_