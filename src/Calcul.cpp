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

mat Calcul::calculPolynomeHermite()
{
    mat H(n_max, z.n_elem);

    for (int i = 0; i < z.n_elem; ++i)
    {
        H(0, i) = 1;
    }
    rowvec h2 = rowvec(z.n_elem);
    h2 = 2 * z;

    H.row(1) = h2;

    for (int i = 2; i < n_max; i++)
    {
        rowvec hn = rowvec(z.n_elem);
        hn = h2 % H.row(i-1) - (2*i) * H.row(i-2);
        H.row(i) = hn;
    }

    return H;
}

Calcul::~Calcul()
{

}
