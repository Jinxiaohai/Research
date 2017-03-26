#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>

#include "arrays.cc"
#include "sf.cc"
#include "random.cc"
#include "misc.cc"
#include "parametermap.cc"

#define CHECKMULTDIV
#define CHECKXCONV

using namespace std;

int main(){
  const double PI=4.0*atan(1.0);
  CCHCalc chcalc;
  CRandom random(-1234);
  int lx,ly,lz,m,L,LMAX=6,i,ntest=10;
  complex<double> AfromY;
  complex<double> ci(0.0,1.0);
  double A;
  double theta,phi,x,y,z;
  
  lx=0; ly=3; lz=0;

  for(i=0;i<ntest;i++){
    theta=acos(1.0-2.0*random.ran());
    phi=2.0*PI*random.ran();

    A=chcalc.GetAFromThetaPhi(lx,ly,lz,theta,phi);
    AfromY=0.0;

    AfromY-=ci*sqrt(PI/35)*SpherHarmonics::Ylm(3,3,theta,phi);
    //AfromY-=sqrt(2*PI/105)*SpherHarmonics::Ylm(3,2,theta,phi);
    AfromY-=ci*(1.0/5.0)*sqrt(3*PI/7)*SpherHarmonics::Ylm(3,1,theta,phi);
    //AfromY-=(2.0/5.0)*sqrt(PI/7)*SpherHarmonics::Ylm(3,0,theta,phi);
    AfromY-=ci*(1.0/5.0)*sqrt(3*PI/7)*SpherHarmonics::Ylm(3,-1,theta,phi);
    //AfromY-=sqrt(2*PI/105)*SpherHarmonics::Ylm(3,-2,theta,phi);
    AfromY-=ci*sqrt(PI/35)*SpherHarmonics::Ylm(3,-3,theta,phi);

    printf("theta=%5.3f phi=%5.3f : A=%8.5f, AfromY=(%8.5f,%8.5f), error=%8.5f\n",
	   theta,phi,
	   A,real(AfromY),imag(AfromY),
	   sqrt(pow(real(AfromY)-A,2)+pow(imag(AfromY),2)));
  }

  return 0;
}
