#include <iostream>
#include <armadillo>
#include "Calcul.h"
#include "DataInterface.h"
#include "Miscellaneous.h"
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

int main()
{
  //orthonormalityTest();
  //HermiteTest();
  //MatrixTest();
  //WnTest();
  //energyTest();
  writingSols();
  return 0;
}
