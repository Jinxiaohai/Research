#ifndef INCLUDE__SF_H
#define INCLUDE__SF_H
#include <gsl/gsl_sf.h>
#include <cstdlib>
#include <cmath>
#include <complex>
#include "misc.h"

using namespace std;

namespace SpherHarmonics{
  double legendre(int ell,double ctheta);
  complex<double> Ylm(int ell,int m,double theta,double phi);
};

namespace Bessel{
  double J0(double x);
  double J1(double x);
  double Jn(int n, double x);

  double K0(double x);
  double K1(double x);
  double Kn(int n, double x);

  double Y0(double x);
  double Y1(double x);
  double Yn(int n, double x);

  double I0(double x);
  double I1(double x);
  double In(int n, double x);

  double j0(double x);
  double j1(double x);
  double jn(int n, double x);

  double y0(double x);
  double y1(double x);
  double yn(int n, double x);

  complex<double> h0(double x);
  complex<double> h1(double x);
  complex<double> hn(int n, double x);

  complex<double> hstar0(double x);
  complex<double> hstar1(double x);
  complex<double> hstarn(int n, double x);

};

namespace CoulWave{

  void GetFG(int L,double x,double eta,double *FL,double *GL);
  void GetFGprime(int L,double x,double eta,double *FL,double *GL,
		  double *FLprime,double *GLprime);
  complex<double> CWincoming(int ell,double x,double eta);
  complex<double> CWoutgoing(int ell,double x,double eta);
  complex<double> cgamma(complex<double> cx);
  void phaseshift_CoulombCorrect(int ell,double q,double eta,
				 double &delta,double &ddeltadq);

};

class CCHCalc{
 public:
  CCHCalc();
  ~CCHCalc();
  // Cartesian Harmonics
  double GetAFromE(int lx,int ly,int lz,double ex,double ey,double ez);
  double GetAFromThetaPhi(int lx,int ly,int lz,double theta,double phi);

  // Returns Overlap \int d\Omega/4pi A_(ell) A(ell')
  double GetOverlap(int lx,int ly,int lz,int lxprime,int lyprime,int lzprime);
  // Same, only does not use pre-stored values, so it is slower if you 
  // calculate many overlaps
  double GetOverlap0(int lx,int ly,int lz,int lxprime,int lyprime,int lzprime);

  // e_x^l_x * e_y^l_y * e_z^l_z
  double GetMFromE(int lx,int ly,int lz,double ex,double ey,double ez);
  double GetMFromThetaPhi(int lx,int ly,int lz,double theta,double phi);

  // Some commonly used functions that are sped up by using stored arrays
  double Factorial(int n);
  double DoubleFactorial(int n);
  double Binomial(int lx,int ly);
  double Trinomial(int lx,int ly,int lz);

 private:
  static int LMAXFACT; // Calculations will bomb if L>LMAXFACT
  static double *fact;
  static double *doublefact;
  static double **binomial;
  static double *****overlap;
  static int INITIALIZED;
  void InitStaticData();
  void ClearStaticData();

  void iswitch(int &i,int &j);
  void overlapinit(int lx,int ly,int lz);

};

#endif

