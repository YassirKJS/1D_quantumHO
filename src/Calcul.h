/**
 * @file Calcul.h
 *
 * Header of Calcul class
 *
 */

#ifndef DEF_CALCUL
#define DEF_CALCUL
#include <armadillo>
#include "Miscellaneous.h"
using namespace arma;

/**
 * @class Calcul
 *
 * This Class computes the Hermite Polynomial and the Wave function
 *
 */

class Calcul
{
    public:
        /// Default constructor of the object
        Calcul();
        /// Constructor of the object that initializes both the n_max and the vector z
        Calcul(int, rowvec);

        mat calculWn();
        /// This function computes the Hermite Polynomial using this recurrence relation
        mat calculPolynomeHermite();

        /// The getter of the n_max property
        int getN();
        /// The setter of the n_max property
        void setN(int);
        /// The getter of the Z property
        void setZ(rowvec z);
        /// The setter of the Z property
        rowvec getZ();

        /// Destructor of the object
        ~Calcul();

        mat calculEnergy(int n, double a, double b, int N);
    private:
        int n_max; /// n_max is the maximum value of n of which we have to compute \f[\psi\f]
        rowvec z; /// z is the entry vector for which we have to compute the Hermite polynomial and the wave function
};

#endif