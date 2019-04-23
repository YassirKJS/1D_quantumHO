#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include "Orthonormality.h"
#include "Calcul.h"
#include "Miscellaneous.h"


/** Checking the othonormality of the wave solutions
*
* This function checks the orthonormality of the wave solutions for two parameters i and j given
*
* @param n_max
* @param j
* @param k
*
*
* @return the result of the orthonormality of the wave solutions
*/

long double Orthonormality::quad(int n_max, int j, int k)
{
	long double sum = 0.0;
	arma::vec Gh;
	arma::rowvec gh;

	arma::mat Hnm;
	arma::mat res;

	arma::rowvec gp =
	{
		-5.38748089, -4.60368245, -3.94476404, -3.34785457, -2.78880606,
		    -2.254974, -1.73853771, -1.23407622, -0.73747373, -0.24534071,
		    0.24534071,  0.73747373,  1.23407622,  1.73853771,  2.254974,
		    2.78880606,  3.34785457,  3.94476404,  4.60368245,  5.38748089
	    };

	arma::mat wgp = {{
			2.22939365e-13, 4.39934099e-10, 1.08606937e-07, 7.80255648e-06,
			2.28338636e-04, 3.24377334e-03, 2.48105209e-02, 1.09017206e-01,
			2.86675505e-01, 4.62243670e-01, 4.62243670e-01, 2.86675505e-01,
			1.09017206e-01, 2.48105209e-02, 3.24377334e-03, 2.28338636e-04,
			7.80255648e-06, 1.08606937e-07, 4.39934099e-10, 2.22939365e-13
		}
	};

	Calcul *cal = new Calcul(n_max, gp);
	arma::mat Wn = cal->calculWn().t();

	for(int n = 0; n < n_max; n++)
	{
		for(int m = n; m < n_max; m++)
		{
			Gh = Wn.col(n) % Wn.col(m);
			gh = Gh.t() / exp(-gp % gp);
			sum = 0.0;
			for(int i = 0; i < 20; i++)
			{
				sum += wgp(i) * gh(i);
			}
			if(n == j && m == k)
			{
				return sum;
				break;
			}
		}
	}
	return 0.0;
}
