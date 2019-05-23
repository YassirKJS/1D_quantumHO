/*! \mainpage <center style="color: #4182ea;">1D Quantum Harmonic Oscillator</center>
 *
 * In this project, I computed the solutions of the <span style="color: #4182ea;">1D Quantum Harmonic Oscillator</span> and I will check some of their properties.<br />
 * <strong>Schrödinger Equation</strong> :
 * <div style="font-size: 30px;"><center>\f$\hat{H}_{(z)}\psi_n(z) = E_n\psi_n(z)\f$</center></div>
 * with the 1D-Hamiltonian and 1D-momentum operators defined as <br/>
 * <div><center>\f$\hat{H}_{(z)}\equiv \frac{\hat{p}_{(z)}^2}{2m} + \frac{1}{2}m\omega^2\hat{z}^2, \hspace{5mm}\hat{p}_{(z)}\equiv -i\hbar\frac{\partial}{\partial z}.\f$</center></div>
 *
 * <div>Here are the solutions of the Schrödinger equation that I have to compute</div>
 * <center>\f$\psi_n(z) = \frac{1}{\sqrt{2^n n!}}\left(\frac{m\omega}{\pi\hbar}\right)^{1/4}e^{-\frac{m\omega z^2}{2\hbar}}H_n\left(\sqrt{\frac{m\omega}{\hbar}} . z\right).\f$</center>
 * The Hermmite polynomial \f$ H_n \f$ will be calculated using the recurence relation below :
 * <div>
 *  <center>\f$H_0(z) = 1\f$</center>
 * <center>\f$H_1(z) = 2z\f$</center>
 * <center>\f$\forall n\ge 1, H_{n+1}(z) = 2zH_n(z)-2nH_{n-1}(z).\f$</center>
 * </div>
 *
 * So after I computed the solutions I plotted the results using R.
 * <img src= "Rplot.png" width="800" />
 * <h3> Orthonormalité </h3>
 * <div style="margin: 15px;">
 * In order  to test our results I will compute the orthonomality. <br/>
 * <center>\f$\forall (m,n), \int \psi^*_m(z)\psi_n(z) dz = \delta_{mn}.\f$</center><br/>
 * I will use the quadrature rule of <strong>Carl Friedrich Gauss</strong> that allows us to approximate the definite integral of a function<br />
 * <center>\f$\int_{-1}^1 f(x) dx \simeq \sum_{i=0}^{n-1}w_if(x_i)\f$</center>
 * There are many quadrature but in our case I will use the <strong>Gauss-Hermite</strong> quadrature which looks as follow : <br />
 * <center><table style="background-color: rgb(0, 136, 192); border-collapse: collapse;" width="50%" border="1">
 *   <tr>
 *     <th>a</th>
 *     <th>b</th>
 *     <th>\f$w(x)\f$</th>
 *      <th>Quadrature</th>
 *     <th>Associated polynomial</th>
 *   </tr>
 *   <tr>
 *      <td style="background-color:white;" align="center">\f$-\infty\f$</td>
 *      <td style="background-color:white;" align="center">\f$+\infty\f$</td>
 *      <td style="background-color:white;" align="center">\f$\(e^{-x^2}\)\f$</td>
 *      <td style="background-color:white;" align="center"><strong>Gauss-Hermite</strong></td>
 *      <td style="background-color:white;" align="center">Hermite</td>
 *    </tr>
 * </table></center>
 * </div>
 */
#include <iostream>
#include <armadillo>
#include "Calcul.h"
#include "DataInterface.h"
#include "Orthonormality.h"
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
