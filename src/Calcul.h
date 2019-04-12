/**
 * @file Calcul.h
   @
 */
#ifndef DEF_CALCUL
#define DEF_CALCUL

#include <armadillo>
using namespace arma;

//! Calcul Class
/*!
    An instance of this class contains the n_max and a vector z
    This Class's main purpose is to compute the solution of the 1D-HO Schrödinger equation
*/
class Calcul
{
public:
    Calcul(); /// Default constructor of the object
    Calcul(int,rowvec); /// Constructor of the object that initializes both the n_max and the vector z

    mat calculWn();
    mat calculPolynomeHermite(); /// This function computes the Hermite Polynomial using this recursive relation

    int getN();
    void setN(int);


    ~Calcul();
private:
    int n_max;
    rowvec z;
};

#endif
