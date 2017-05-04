#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#define double_complex complex
#define pi 3.14159265358979323844
#define EULER 0.5772156649015328606
double_complex ci(0.0,1.0);

#define STRONG_RMAX 4.0
#define COULOMB
#define STRONG_NRMAX 100
#define MASS1 938.3
#define MASS2 938.3
#define Q1Q2 1
#define POTENTIAL vreid

static int STRONG_L[1]={0};

double_complex cw2_big_r(int l,double r,double eta);
double_complex cw2_small_r(int l,double r,double eta);
double_complex old_cw2_small_r(int l,double r,double eta);
double_complex cgamma(double_complex c);
double dgamma(int mm);
double vreid(double r,int ipart);

void main(){
  double delr,r00,rmax,mom,phase,phase0,phase1,phase2;
  double vv,r1,r2,r0,vcs;
  double b,eta;
  const double rmass=MASS1*MASS2/(MASS1+MASS2);
  double_complex  cf[STRONG_NRMAX],cf0[STRONG_NRMAX],cfin[2*STRONG_NRMAX+1];
  double_complex cy0,cy1,cy2,ceip,ctest,ctest0,ctest1,ctest2,oldctest;
  int kk,j,mm,nmax,istrong,l,ipart;
  //Wave functions are found as a function of kr, referred to as
  //r here.  Values are calculated for delr=.d0005*krmax,
  //krmax=40fm*mom
  //other value is passed on to the main prog.  
  //So the main prog. uses
  //a mesh of delr=.d001*krmax
  nmax=2*STRONG_NRMAX;
  ipart=0;
  l=STRONG_L[ipart];

  printf("what is mom?\n");
  scanf("%lf",&mom);
  rmax=STRONG_RMAX*mom/197.323;
  delr=-rmax/(double)nmax;
  r00=rmax;
  //-----------------
  //if(r00<6&&eta&&1.0E-5){
  //r00=-delr*int(-6.0/delr)
  //}
  //-----------------
#ifdef COULOMB
  eta=(double)Q1Q2*(rmass/137.036)/mom;
  printf("eta=%g\n",eta);
#else
  eta=0.0;
#endif
  //+++++++++++++++ W & WO POTENTIAL LOOP +++++++++++++++++
  //We perform calc with and with out the Reid pot.
  //We are only interested in the correction to the wave function
  //due to interactions. The first
  //time through we turn off the Reid pot.(istrong=1,off) (istrong=2,on)
  for(istrong=0;istrong<2;istrong++){
    r0=r00;
    r1=r0+delr;
    //This initiates an incoming coul. wave function at
    //two points at r=10.  See Messiah's appendix on incoming partial
    //Coul. waves for explanation.  We then calc. inward until r=0.
#if ! defined COULOMB
    cy0=cw2_big_r(l,r0,eta)*exp(-ci*r0);
    cy1=cw2_big_r(l,r1,eta)*exp(-ci*r1);
#else
    if(r0>10.0){
      printf("do the big one\n");
      cy0=cw2_big_r(l,r0,eta)*exp(-ci*r0);
      cy1=cw2_big_r(l,r1,eta)*exp(-ci*r1);
    }
    else{
      printf("do the small one\n");
      cy0=cw2_small_r(l,r0,eta);
      cy1=cw2_small_r(l,r1,eta);
    }
#endif
    
    cy0=cy0*pow(ci,l);
    cy1=cy1*pow(ci,l);
    printf("cy1=%g %g\n",real(cy1),imag(cy1));


    cfin[nmax]=cy0;
    cfin[nmax-1]=cy1;
    //Given the wave func. at two pts the wave func. can be found at
    //third pt. from diff. eq.
    for(j=nmax-2;j>0;j=j-1){
      r2=r1+delr;
      vcs=1.0;
      if(istrong==0){
	vv=0.0;
      }
      else{
	vv=(2.0*rmass/(mom*mom))*POTENTIAL(197.323*r1/mom,ipart);
      }
      /* For screened coulomb */
      //rvcs=(r1/mom)*197.323;
      //vcs=(rvcs/4.d0)**3
      b=1.0-((double)(l*(l+1))/pow(r1,2))-(2.0*eta*vcs/r1)-vv;
      cy2=(-b*cy1*pow(delr,2)+2*cy1-cy0);
      
      /*ctest0=cw2_small_r(l,r0,eta);
	ctest1=cw2_small_r(l,r1,eta);
	ctest2=cw2_small_r(l,r2,eta);
	ctest=(ctest2-2.0*ctest1+ctest0)+b*ctest1*delr*delr;
	ctest=ctest/(ctest2-ctest0);
	printf("r=%6g ------ cy1=%6g,%6g  ctest=%6g %6g\n",
	r1*197.323/mom,real(cy1),imag(cy1),real(ctest),imag(ctest));
	printf("r2=%4g psi2=%6g,%6g should be:%6g,%6g ratio=%6g,%6g\n",
	r2*197.323/mom,real(cy2),imag(cy2),real(ctest2),imag(ctest2),
	real(cy2/ctest2),imag(cy2/ctest2)); */
      
      cfin[j]=cy2;
      r0=r1;
      r1=r2;
      cy0=cy1;
      cy1=cy2;
    }
    //phase0 is the phase of the wavefunction at r=0.  
    //Subtracting the double_complex
    //conjugate of the wavefunction *exp(-2i*phase0) 
    //from the wavefunction
    //yields a sol. which satisfies the bound. conditions.Multiply by i
    //to correspond to the correct incoming partia-wave phase(
    //See Messiah).
    //finally, when coming through after doing it for both with and 
    //withoutinteractions, take the difference.
    phase1=real(-ci*log(cfin[1]/abs(cfin[1])));
    phase2=real(-ci*log(cfin[2]/abs(cfin[2])));
    phase0=2.0*phase1-phase2;
    ceip=exp(ci*phase0);
    //Note that we send back to the main prog. the answer with a less
    //dense mesh
    j=0;
    for(mm=1;mm<2*STRONG_NRMAX;mm=mm+2){
      if(istrong==0){
	cf0[j]=ci*(cfin[mm]-conj(cfin[mm])*pow(ceip,2))/2;
      }
      if(istrong==1){
	cf[j]=ci*(cfin[mm]-conj(cfin[mm])*pow(ceip,2))/2;
	r0=-delr*((double)mm);
      }
      j=j+1;
    }
    if(istrong==0) phase=phase0;
    if(istrong==1) {
      printf("%6g  %g %g\n",mom,phase*180.0/pi,(phase0-phase)*180.0/pi);
    }
  }
  for(j=0;j<STRONG_NRMAX;j++){
    r1=-((double)j+0.5)*2.0*delr*197.323/mom;
    printf("%5g   %g,%g\n",r1,real(cf[j]-cf0[j])/r1,imag(cf[j]-cf0[j])/r1);
  }
}
//**********************************
double_complex cw2_big_r(int l,double r,double eta){
  double_complex z,a,b,f1,top1,top2,bot,delf,answer;
  double arg;
  int n,j;
  z=-2*r*ci;
  a=(double)l+1.0+ci*eta;
  b=2.0*((double)l+1.0);
  f1=1.0;
  top1=1.0;
  top2=1.0;
  bot=1.0;
  n=5;
  for(j=1;j<=n;j++){
    top1=top1*((double)j-a);
    top2=top2*((double)j+b-a-1.0);
    bot=bot*(double)j;
    delf=(top1*top2)/(bot*pow(z,j));
    f1=f1+delf;
  }
  arg=eta*log(2.0*r);
  arg=arg-2.0*pi*floor(arg/(2.0*pi));
  answer=exp(ci*arg)*f1;
  return answer;
}
//*************************************
double_complex cw2_small_r(int l,double r,double eta){
  //The notation is like Gr. + R, page 1063.
  //The Coulomb wave function is the same as W(i*eta,l+1/2,2*i*rho)
  double_complex factor,lterm,fact1,fact2;
  double_complex psi1,psi2,psi3,cx,sum1,sum2,delsum1,delp,cdcon1,cdcon2,answer;
  int k;
#ifdef COULOMB
  cdcon1=cgamma(-(double)(l)-ci*eta);
  cdcon2=cgamma((double)(l+1)-ci*eta);
#else
  cdcon1=dgamma(-l);
  cdcon2=dgamma(l+1);
#endif
  factor=pow(-1,2*l+1)*pow(2*ci*r,l+1)*exp(-ci*r)/(cdcon1*cdcon2);
  psi1=-EULER;
  psi2=-EULER;
  for(k=1;k<=2*l+1;k++){
    psi2=psi2+1.0/(double)k;
  }
  cx=(double)(l+1)-ci*eta;
  psi3=-EULER-(1.0/cx)+cx*(pi*pi/6.0);
  for(k=1;k<100000;k++){
    delp=-cx*cx/((double)(k*k)*(cx+(double)k));
    psi3=psi3+delp;
    if(abs(delp)<1.0E-12) goto CONVERGE1;
  }
  printf("never escaped loop1 in cw2_small_r!\n");
CONVERGE1:
  lterm=log(2*ci*r);
  fact1=cdcon2/dgamma(2*l+2);
  sum1=fact1*(psi1+psi2-psi3-lterm);
  for(k=1;k<=10000;k++){
    fact1=fact1*(2*ci*r)*((double)(l+k)-ci*eta)/((double)k*(double)(2*l+1+k));
    psi1=psi1+1.0/(double)k;
    psi2=psi2+1.0/(double)(2*l+1+k);
    psi3=psi3+1.0/((double)(k-1)+cx);
    delsum1=fact1*(psi1+psi2-psi3-lterm);
    sum1=sum1+delsum1;
    if(abs(delsum1)<1.0E-15) goto CONVERGE2;
  }
  printf("never escaped loop2 in cw2_small_r!\n");
CONVERGE2:
  fact2=dgamma(2*l+1)*cdcon1/pow(-2*ci*r,2*l+1);
  sum2=fact2;
  for(k=1;k<=2*l;k++){
    fact2=fact2*((double)(k-l-1)-ci*eta)*
      (-2.0*ci*r)/((double)(k)*(double)(2*l-k+1));
    sum2=sum2+fact2;
  }
  sum1=factor*sum1;
  sum2=factor*sum2;
  answer=(sum1+sum2)*exp(pi*eta/2.0);
  return answer;
}
//******************************************
double dgamma(int mm){
  //This calc.s gamma functions which are in the form gamma(n)
  //where n is an int > 0.
  double cg;
  int j;
  cg=1.0;
  if(mm<1) {
    for(j=1;j<=-mm+1;j++){
      cg=cg/(1.0+(double)(-j));
    }
  }
  if(mm>1){
    for(j=1;j<=mm-1;j++){
      cg=cg*((double)j);
    }
  }
  return cg;
}
//******************************************
double_complex cgamma(double_complex c){
  //This calc.s gamma functions which are in the form gamma(n+i*y)
  //where n is an int and y is real.
  double_complex cg,cphase;
  int mm,j;
  double x,y,phase,delp,cgmag;
  x=real(c);
  y=imag(c);
  phase=-EULER*y;
  for(j=1;j<=100000;j++){
    delp=(y/(double)j)-atan(y/(double)j);
    phase=phase+delp;
    if(fabs(delp)<1E-10) goto CGAMMA_ESCAPE;
  }
  printf("oops not accurate enough, increase jmax\n");
CGAMMA_ESCAPE:
  phase=phase-2.0*pi*floor(phase/(2.0*pi));
  cphase=exp(ci*phase);
  cgmag=sqrt(pi*y/sinh(pi*y));
  mm=(int)floor(x+0.5);
  cg=cgmag*cphase;
  if(mm<1){
    for(j=1;j<=-mm+1;j++){
      cg=cg/(1.0+(double)(-j)+ci*y);
    }
  }
  if(mm>1) {
    for(j=1;j<=mm-1;j++){
      cg=cg*((double)(j)+ci*y);
    }
  }
  return cg;
}
double vreid(double r,int ipart){
  double pmux,f1,f4,f7,f2,f3,vr;
  if(ipart==0){
    /* l=0 */
    pmux=r*0.7; /*0.7 is 1/mpi in fm */
    f1=exp(-pmux);
    f4=pow(f1,4);
    f7=f4*pow(f1,3);
    vr=-10.463*f1/pmux-1650.6*f4/pmux;
    vr=vr+6484.2*f7/pmux;
  }
  if(ipart>0){
    /* l=1, j=0,1,2 */
    pmux=r*0.7; /*0.7 is 1/mpi in fm */
    f1=exp(-pmux);
    f2=f1*f1;
    f4=f2*f2;
    f3=f2*f1;
    vr=(1.0+2.0/pmux+2.0/(pmux*pmux))*f1;
    vr=vr-(8./pmux+2./(pmux*pmux))*f4;
    vr=10.463*vr/pmux;
    vr=vr-135.25*f2/pmux+472.81*f3/pmux;
  }
  return vr;
}
