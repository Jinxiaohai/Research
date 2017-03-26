#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <complex>
using namespace std;

#include "misc.cc"
#include "sf.cc"

using namespace std;

int main(){
  int n,ell,m;
  double theta,phi,ctheta;
  complex<double> cans;
  double ans;
  
 TRYNEW:
  printf("What are theta and phi? : ");
  scanf("%lf %lf",&theta,&phi);
  ell=3; m=-1;
  ctheta=cos(theta);

  cans=SpherHarmonics::Ylm(ell,m,theta,phi);
  ans=SpherHarmonics::legendre(ell,ctheta);
  printf("Y(l=%d,m=%d,theta=%g)=(%g,%g), P(l=%d,ctheta=%g)=%g\n",
	 ell,m,theta,real(cans),imag(cans),ell,ctheta,ans);

  printf( "Enter 0 to quit, other to continue : ");
  scanf("%d",&n);
  if(n!=0) goto TRYNEW;

  return 0;
  
}
