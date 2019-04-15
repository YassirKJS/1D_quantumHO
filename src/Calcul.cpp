#include "Calcul.h"
#import <math.h>

Calcul::Calcul()
{

}

Calcul::Calcul(int n_max, rowvec z)
{
    this->n_max = n_max;
    this->z = z;
}

int Calcul::getN()
{
    return this->n_max;
}
void Calcul::setN(int n_max)
{
    this->n_max = n_max;
}
rowvec Calcul::getZ()
{
    return this->z;
}
void Calcul::setZ(rowvec z)
{
    this->z = z;
}

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

int factorial(int n)
{
    if (n == 0)
        return 1;
    else if (n == 1)
        return 1;
    else if (n==2)
        return 2;
    else
        return n * factorial(n-1); // recursive call to factorial()
}

mat Calcul::calculWn()
{
    rowvec exp_term = exp(-(z%z)/2);
    mat H = this->calculPolynomeHermite();
    mat W(n_max, z.n_elem);
    for (int i = 0; i < this->n_max; ++i)
    {
        double first_term = (1/sqrt(pow(2, i)*factorial(i))) * pow(1/datum::pi, 0.25);
        W.row(i) = first_term * (exp_term % H.row(i));
    }
    return W;
}

Calcul::~Calcul()
{

}
