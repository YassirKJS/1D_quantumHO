/*
 * @file Calcul.cpp
 *
 * Implementation of the class Calcul
 *
 */
#include <stdlib.h>
#include "Calcul.h"
#include <math.h>
#include <cmath>


/** Default Constructor of Calcul class
 *
 *
 */
Calcul::Calcul()
{

}

/** Constructor of Calcul class with given values
 *
 * @param n_max is the maximum value of n of which we have to compute \f[\psi\f]
 * @param z is the entry vector
 */
Calcul::Calcul(int n_max, rowvec z)
{
    this->n_max = n_max;
    this->z = z;
}

/** Getter of the n_max property
 *
 * @return This function returns the n_max value
 *
 */
int Calcul::getN()
{
    return this->n_max;
}
/** Setter of the n_max property
 *
 * This function sets the value of n_max
 *
 * @param n_max
 *
 */
void Calcul::setN(int n_max)
{
    this->n_max = n_max;
}

/** Getter of the z property
 *
 * @return This function returns the current vector z;
 *
 */
rowvec Calcul::getZ()
{
    return this->z;
}

/** Setter of the z property
 *
 * This function sets the value of n_max
 *
 * @param z
 *
 */
void Calcul::setZ(rowvec z)
{
    this->z = z;
}

/** Computing the Hermite Polynomial
 *
 *
 * This function computes the Hermite Polynomial using the reccurence formula
 *
 * \#Definition(physicists version):
 * \f[\forall n\ge 0, H_n(z)\equiv (-1)^n e^{z^2} \frac{d^n}{dz^n}\left( e^{-z^2} \right).\f]
 *
 * \#Reccurence relation:
 * \f[H_0(z) = 1\f]
 * \f[H_1(z) = 2z\f]
 * \f[\forall n\ge 1, H_{n+1}(z) = 2zH_n(z)-2nH_{n-1}(z).\f]
 *
 * @return a matrix containing the values of the function
 */


mat Calcul::calculPolynomeHermite()
{
    mat H(n_max, z.n_elem);

    if(n_max == 0)
    {
        H = z.ones(size(z));
    }

    else
    {
        if(n_max == 1)
        {
            return z.for_each([](arma::mat::elem_type & val)
            {
                val = 2 * val;
            });
        }
        else
        {
            for(int i = 0; i < z.n_elem; ++i)
            {
                H(0, i) = 1;
            }
            rowvec h2 = rowvec(z.n_elem);
            h2 = 2 * z;

            H.row(1) = h2;

            for(int i = 2; i < n_max; i++)
            {
                rowvec hn = rowvec(z.n_elem);
                hn = h2 % H.row(i - 1) - (2 * i) * H.row(i - 2);
                H.row(i) = hn;
            }
        }
    }
    return H;
}


/** Computing the solution of the wave function
 *
 *
 * This function computes the wave function solutions
 *
 * \#The 1D-HO Schrödinger equation:
 * \f[\left(\frac{\hat{p}_{(z)}^2}{2m} + \frac{1}{2}m\omega^2\hat{z}^2\right)\psi_n = E_n\psi_n.\f]
 *
 * \#The analytic solutions form of the wave function:
 * \f[\psi_n(z) = \frac{1}{\sqrt{2^n n!}}\left(\frac{m\omega}{\pi\hbar}\right)^{1/4}e^{-\frac{m\omega z^2}{2\hbar}}H_n\left(\sqrt{\frac{m\omega}{\hbar}} . z\right).\f]
 *
 * @return a matrix containing the values for each
 */
mat Calcul::calculWn()
{

    Miscellaneous misc;
    rowvec exp_term = exp(-(z % z) / 2);
    mat H = this->calculPolynomeHermite();
    mat W(n_max, z.n_elem);
    for(int i = 0; i < this->n_max; ++i)
    {
        double first_term = (1 / sqrt(pow(2, i) * misc.factorial(i))) * pow(1 / datum::pi, 0.25);
        W.row(i) = first_term * (exp_term % H.row(i));
    }
    return W;

}

/**
*
* Destructor of the object Calcul
*/

Calcul::~Calcul()
{

}

mat Calcul::calculEnergy(int n, double a, double b, int N)
{
  mat E;
  Miscellaneous misc;
  double pas = (b - a) / N;
  rowvec Z = arma::linspace<rowvec>(a, b + 2 * pas, N + 2);
  this->z = Z;
  this->n_max = n;
  mat Psi = this->calculWn().t();
  mat ZZ(n, N + 2);
  for(int i = 0; i < ZZ.n_rows; i++)
  {
    ZZ.row(i) = z;
  }
  ZZ = ZZ.t();

  mat dPsi = misc.deriv(Psi.rows(0, N), Psi.rows(1, N + 1), ZZ.rows(0, N), ZZ.rows(1, N + 1));
  mat d2Psi = misc.deriv(dPsi.rows(0, N - 1), dPsi.rows(1, N), ZZ.rows(0, N - 1), ZZ.rows(1, N));
  mat ZN = ZZ.rows(0, N - 1);
  mat PsiN = Psi.rows(0, N - 1);

  //E = -1 * h * h * d2Psi / (2 * m * PsiN) + m * w * w * ZN % ZN / 2;
  E = -1 * 1 * 1 * d2Psi / (2 * 1 * PsiN) + 1 * 1 * 1 * ZN % ZN / 2;
  return E;
}