#include <iostream>
#include <armadillo>
#include "Calcul.h"
#include "DataInterface.h"
#include "Miscellaneous.h"
#include "Orthonormality.h"
using namespace arma;
using namespace std;

void test()
{
  //mat z = rowvec({-2, -1, 0, 1, 2});
  Calcul *c = new Calcul();
  DataInterface d;
  int n = 5;
  rowvec z = d.read(n);

  z.print();
  cout << n << endl;

  c->setN(n);
  c->setZ(z);

  d.write(z, c->calculWn());
  // d.write(z, c->calculPolynomeHermite());
}

void HermiteTest()
{
  cout << "starting HermiteTest()" << endl;
  mat z = rowvec({-2, -1, 0, 1, 2});
  Calcul *c1 = new Calcul(3, z);
  Calcul *c2 = new Calcul(4, z);
  mat H1 = c1->calculPolynomeHermite();
  mat H2 = c2->calculPolynomeHermite();
  cout << "H(3,z):" << endl << H1 << endl;
  cout << "H(4,z):" << endl << H2 << endl;
}

void MatrixTest()
{
  mat A = { {1, 2, 3},
    {4, 5, 6},
    {7, 8, 9}
  };
  A.print();

  A.row(0).print();
  A.row(1).print();
}

void WnTest()
{
  Calcul *c = new Calcul();
  int n = 1;
  mat z = rowvec({-2, -1, 0, 1, 2});
  z.print();
  cout << n << endl;

  c->setN(n);
  c->setZ(z);

  mat Wn = c->calculWn();
  Wn.print();
}

void writingSols()
{
  //test();
  mat z = rowvec({-2, -1, 0, 1, 2});
  Calcul *c = new Calcul(3, z);
  mat W = c->calculWn();
  DataInterface d;
  d.write(z, W);
}

void orthonormalityTest()
{
  std::cout << "Testing Orthonormality n = " << 4 << " m = " << 4 << " quad: " << Orthonormality::quad(6, 1, 1) << std::endl;
}

void energyTest()
{
  mat E;
  int i;
  int N = 1000000;
  int n = 5;

  // Miscellaneous misc;
  double b = 2;
  double a = -2;

  Calcul *cal = new Calcul();

  for(i = 2; i <= n; i++)
  {
    cal->setN(i);
    mat res;
    res = cal->calculEnergy(i, a, b, N);
  }

  mat E2 = cal->calculEnergy(2, -2, 2, 1000000);
  cout << "E2: " << E2.rows(2, 6) << endl;
  mat m1 = {{0.5, 0.5, 0.5, 0.5, 0.5}, {1.5, 1.5, 1.5, 1.5, 1.5}};
}

int main()
{
  energyTest();
  //orthonormalityTest();
  //HermiteTest();
  //MatrixTest();
  //WnTest();
  //writingSols();
  return 0;
}
