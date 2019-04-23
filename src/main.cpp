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

