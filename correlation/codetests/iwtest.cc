#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
using namespace std;
#include "coral.cc"

int main(){
  const double PI=3.14159265358979323844;
  const double HBARC=197.3269602,ROOT2=sqrt(2.0);
  double delta=0.23456;
  int ell=1,ellmax=1,q1q2;
  complex<double> psi,hstar,h,hstarq0,hq0,ci(0.0,1.0);
  double r,x,q=34.0,epsilon=1.0,eta,e1,e2;
  double r0=1.0,rf=4.0,delr=0.001,iw=0.0,iwq0=0.0;
  CPartWave *partwave,*partwaveq0;
  CWaveFunction *wf;
  wf=new CWaveFunction();
  printf("Enter q and delta and ell\n");
  scanf("%lf %lf %d",&q,&delta,&ell);
  //q=35.0;
  //delta=0.4;
  //ell=0;
  q1q2=-1;
  e1=e2=sqrt(q*q+139.57*139.57);
  eta=double(q1q2)*e1*e2/((e1+e2)*137.036*q);
  //eta=1.0E-8;
  printf("eta=%g\n",eta);
  partwave=new CPartWave(eta,q1q2,q,ell,epsilon);
  partwaveq0=new CPartWave(0,0,q,ell,epsilon);
  for(r=r0+0.5*delr;r<rf;r+=delr){
    x=q*r/HBARC;
    hstar=partwave->GetPhiIncoming(r)/x;
    h=partwave->GetPhiOutgoing(r)/x;
    hstarq0=partwaveq0->GetPhiIncoming(r)/x;
    hq0=partwaveq0->GetPhiOutgoing(r)/x;
    /*
      printf("h=(%g,%g), hq0=(%g,%g), e^(ix)/x=(%g,%g)\n",
      real(h),imag(h),real(hq0),imag(hq0),
      real(-ci*ceiphi(x)/x),imag(-ci*ceiphi(x)/x));
      printf("hstar=(%g,%g), hstarq0=(%g,%g)\n",
      real(hstar),imag(hstar),real(hstarq0),imag(hstarq0));
    */
    psi=hstar*Misc::ceiphi(-2.0*delta)+h;
    iw+=0.5*x*x*real(psi*conj(psi))*delr/HBARC;
    psi=hstarq0*Misc::ceiphi(-2.0*delta)+hq0;
    iwq0+=0.5*x*x*real(psi*conj(psi))*delr/HBARC;
  }
  printf("IW=%g =? %g, IWq0=%g =? %g\n",iw,
	 wf->GetIW(ell,r0,q,q1q2,eta,delta)-wf->GetIW(ell,rf,q,q1q2,eta,delta),
	 iwq0,
	 wf->GetIW(ell,r0,q,0,0.0,delta)-wf->GetIW(ell,rf,q,0,0.0,delta));
  delete partwave;

  return 0;
  
}
