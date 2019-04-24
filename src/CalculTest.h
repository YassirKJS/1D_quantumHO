/**
 * @file CalculTest.h
 * the Calcul Unit Test
 * In this Test class we will test all the computation methods
 */
#ifndef DEF_CALCUL_TEST
#define DEF_CALCUL_TEST

#include <cxxtest/TestSuite.h>
#include <armadillo>
#include "Calcul.h"
#include "Miscellaneous.h"
#include "Orthonormality.h"

using namespace arma;

class CalculTest : public CxxTest::TestSuite
{
	public:

		/**
		*    Test of the Computation of the Hermite Polynomial
		*/
		void testPolynomeHermite(void)
		{
			mat z = rowvec({-2, -1, 0, 1, 2});
			Calcul *c = new Calcul(3, z);
			mat H = c->calculPolynomeHermite();
			mat expected_H = {{1, 1, 1, 1, 1},
				{-4, -2, 0, 2, 4},
				{12, 0, -4, 0, 12}
			};

			TS_TRACE("********** STARTING THE COMPARISON BETWEEN THE RESULT AND THE EXPECTED ONE ************");

			TS_TRACE("START TEST FOR H0");
			double difference = norm(H.row(0) - expected_H.row(0), 2);
			TS_ASSERT_EQUALS(difference, 0);
			TS_TRACE("END TEST FOR H0");

			TS_TRACE("START TEST FOR H1");
			difference = norm(H.row(1) - expected_H.row(1), 2);
			TS_ASSERT_EQUALS(difference, 0);
			TS_TRACE("END TEST FOR H1");

			TS_TRACE("START TEST FOR H2");
			difference = norm(H.row(2) - expected_H.row(2), 2);
			TS_ASSERT_EQUALS(difference, 0);
			TS_TRACE("END TEST FOR H2");

			TS_TRACE("********** END OF THE COMPARISON BETWEEN THE RESULT AND THE EXPECTED ONE ************");
		}

		/**
		*    Test of the Computation of the wave function
		*/
		void testcalculWn(void)
		{
			mat z = rowvec({-2, -1, 0, 1, 2});
			Calcul *c = new Calcul(3, z);
			mat Wn = c->calculWn();
			mat expected_W = {{0.101654, 0.455581, 0.751126, 0.455581, 0.101654},
				{-0.28752, -0.644288, 0, 0.644288, 0.28752},
				{0.43128, 0, -1.062252, 0, 0.43128}
			};

			TS_TRACE("********** STARTING THE COMPARISON BETWEEN THE RESULT AND THE EXPECTED ONE ************");

			TS_TRACE("START TEST FOR W0");
			double difference = norm(Wn.row(0) - expected_W.row(0), 2);
			TS_ASSERT_DELTA(difference, 0, 0.00001);
			TS_TRACE("END TEST FOR W0");

			TS_TRACE("START TEST FOR W1");
			difference = norm(Wn.row(1) - expected_W.row(1), 2);
			TS_ASSERT_DELTA(difference, 0, 0.00001);
			TS_TRACE("END TEST FOR W1");

			TS_TRACE("START TEST FOR W1");
			difference = norm(Wn.row(2) - expected_W.row(2), 2);
			TS_ASSERT_DELTA(difference, 0, 0.00001);
			TS_TRACE("END TEST FOR W1");

			TS_TRACE("********** END OF THE COMPARISON BETWEEN THE RESULT AND THE EXPECTED ONE ************");
		}

		/**
		*    Testing the factorial method
		*/
		void testFactorial(void)
		{
			Miscellaneous misc;
			TS_TRACE("Testing Factorial");
			int fac = misc.factorial(3);
			TS_ASSERT_EQUALS(fac, 6);
		}

    /**
    *    Testing the Orthonormality
    */

    void testOrthonormality(void)
    {
      TS_TRACE("Testing Orthonormality");
      double q1 = Orthonormality::quad(6, 0, 0);
      double q2 = Orthonormality::quad(6, 1, 1);
      double q3 = Orthonormality::quad(5, 1, 1);
      double q4 = Orthonormality::quad(6, 0, 1);
      double q5 = Orthonormality::quad(6, 2, 3);
      TS_ASSERT_DELTA(q1, 1.0, 0.0000001);
      TS_ASSERT_DELTA(q2, 1.0, 0.0000001);
      TS_ASSERT_DELTA(q3, 1.0, 0.0000001);
      TS_ASSERT_DELTA(q4, 0.0, 0.0000001);
      TS_ASSERT_DELTA(q5, 0.0, 0.0000001);
    }
		
};

#endif
