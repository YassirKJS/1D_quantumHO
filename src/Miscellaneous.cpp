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

arma::mat Miscellaneous::deriv(arma::mat Y1, arma::mat Y2, arma::mat X1, arma::mat X2)
{
  arma::mat derv = (Y2 - Y1) / (X2 - X1);
  return derv;
}
