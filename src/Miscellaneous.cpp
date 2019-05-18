/**
 * @file Miscellaneous.cpp
 *
 * \#Implementation of Miscellaneous functions for calculus
 */

#include "Miscellaneous.h"


/** Computing the factorial of n
*
* This function computes the factorial of n using reccurence
*
* @param n integer
*
* @return the result of the factorial of n
*/
int Miscellaneous::factorial(int n)
{
	if(n == 0)
		return 1;
	else if(n == 1)
		return 1;
	else if(n == 2)
		return 2;
	else
		return n * factorial(n - 1); // recursive call to factorial()
}

/**
 * Computing the derivative of a function
 *
 * this function approximates the fist order derivative
 *
 * @param Y1 vector containing the n values of the image
 * @param Y2 vector containing the n values of the image offset by one value
 * @param X1 vector containing the n points where the image is valued
 * @param X2 vector containing the n points where the image is valued offset by one value
 *
 * @return derv vector containing the values of the derived function
 */
arma::mat Miscellaneous::deriv(arma::mat Y1, arma::mat Y2, arma::mat X1, arma::mat X2)
{
  arma::mat derv = (Y2 - Y1) / (X2 - X1);
  return derv;
}
