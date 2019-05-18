#ifndef DEF_MISCELLANEOUS
#define DEF_MISCELLANEOUS
/**
 * @file Miscellaneous.h
 *
 * Header of the Miscellaneous Class
 *
 */
#include <iostream>
#include <armadillo>

/**
 * @class Miscellaneous
 *
 * \#Miscellaneous functions for calculus
 */
class Miscellaneous
{
	public:
		int factorial(int n);
    arma::mat deriv(arma::mat Y1, arma::mat Y2, arma::mat X1, arma::mat X2);
};

#endif
