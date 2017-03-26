#ifndef __INCLUDE_MISC_H
#define __INCLUDE_MISC_H

#include <complex>
#include <cstdlib>
#include <gsl/gsl_sf.h>
#include "sf.h"
#include "random.h"

using namespace std;

namespace Misc{
  void lorentz(double *u,double *p1,double *p1prime);
  double cgc(double j1,double m1,double j2,double m2,double j,double m);
  bool comparestrings(char *s1,char *s2);
  double triangle(double m0,double m1,double m2);
  void outsidelong(double *pa,double *pb,
		   double &qinv,double &qout,double &qside,double &qlong);
  double GetQinv(double *pa,double *pb);
  double GetRapidity(double *pa);
  double GetDely(double *pa,double *pb);
  
  complex<double> cexp(complex<double> z);
  complex<double> ceiphi(double phi);
  complex<double> cpow(complex<double> z,complex<double> a);

  int iround(double x);

  int cgc_delta (int x, int y);
  double cgc_factorial (double n);
  double cgc_fractorial (double n,double m);
  double oldcgc(double j1,double m1,double j2,double m2,double j,double m);
};

#endif
