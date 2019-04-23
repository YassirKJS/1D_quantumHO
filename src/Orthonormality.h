#ifndef DEF_ORTH
#define DEF_ORTH

#include <armadillo>

/**
 * @class Orthonormality
 *
 * # Orthogonality
 * \f[\forall (m,n), \int \psi^*_m(z)\psi_n(z) dz = \delta_{mn}.\f]
 *
 * # Common weight function
 * \f[\int_{a}^b \omega(x)g(x) dx \simeq \sum_{i=0}^{n-1}w^\omega_ig(x^\omega_i)\f]
 * | \f$a\f$       | \f$b\f$       | \f$\omega(x)\f$                                    | Quadrature      | Associated polynomial |
 * | :-----------: | :-----------: | :------------------------------------------------: | :-------------: | :-------------------: |
 * | \f$−\infty\f$ | \f$+\infty\f$ | \f$e^{-x^2}\f$                                     | Gauss–Hermite   | Hermite               |
 *
 */

using namespace arma;
class Orthonormality
{
	public:
		static long double quad(int n_max, int n, int m);
};

#endif
