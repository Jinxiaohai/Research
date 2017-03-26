#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
using namespace std;
#include <cstdio>
#include <cmath>
#include <gsl/gsl_sf.h>

int main(){
  int L=2,dummy;
  double x,eta=0.01;
  double expF,expG;
  double dphi0F,dphi1F,dphi2F,dphi0G,dphi1G,dphi2G;
  double F0,F1,F2,G0,G1,G2;
  double *fc,*gc;
  fc=new double[L+1];
  gc=new double[L+1];

  printf("eta=%g\n",eta);
  printf("x ____ F(L=0)  G(L=0) ___  F(L=1) G(L=1) ___ F(L=2) G(L=2) ___\n");
  for(x=0.1;x<7;x+=0.05){
    dummy=gsl_sf_coulomb_wave_FG_array(0,L,eta,x,fc,gc,&expF,&expG);
    expF=exp(expF);
    expG=exp(expG);
    F0=expF*fc[0];
    F1=expF*fc[1];
    F2=expF*fc[2];
    G0=expG*gc[0];
    G1=expG*gc[1];
    G2=expG*gc[2];
    printf("%4.2f   %7.3f %7.3f   %7.3f %7.3f   %7.3f %7.3f\n",
	   x,F0,G0,dphi0F,dphi0G,F1,G1,dphi1F,dphi1G,F2,G2,dphi2F,dphi2G);
  }
  delete [] fc;
  delete [] gc;

  return 0;
}








