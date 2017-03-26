#ifndef INCLUDE_SOURCECALC_CC
#define INCLUDE_SOURCECALC_CC
using namespace std;

#include "sourcecalc.h"

CSourceCalc::CSourceCalc(){
}

void CSourceCalc::CalcS(CCHArray *A){
  printf("i'm a dummyCalcS(CCHArray)\n");
  // Dummy function to be over-written by inherited class
}
void CSourceCalc::CalcS(int lx,int ly,int lz,CCHArray *A){
  printf("i'm a dummy CalcS(int,int,int,CCHArray *)\n");
  // Dummy function to be over-written by inherited class
}


void CSourceCalc::ReadSPars(char *sparsfilename){
  parameter::ReadParsFromFile(spars,sparsfilename);
}

void CSourceCalc::NormCheck(CCHArray *A){
  const double PI=4.0*atan(1.0);
  double r,check,DELR;
  int ir,NRMAX;
  NRMAX=A->GetNRADIAL();
  DELR=A->GetRADSTEP();
  check=0;
  for(ir=0;ir<NRMAX;ir++){
    r=(0.5+ir)*DELR;
    check+=r*r*DELR*A->GetElement(0,0,0,ir)*(4.0*PI);
  }
  printf("normalization check = %g\n",check);
}

void CSourceCalc::CalcEffGaussPars(CCHArray *A){
  double xbar,ybar,zbar,x2bar,y2bar,z2bar,r2bar,r,r2,r3,r4,DELR;
  int ir,NRMAX;
  bool XSYM,YSYM,ZSYM;
  NRMAX=A->GetNRADIAL();
  DELR=A->GetRADSTEP();
  XSYM=A->GetXSYM();
  YSYM=A->GetYSYM();
  ZSYM=A->GetZSYM();
  const double PI=4.0*atan(1.0);
  xbar=ybar=zbar=x2bar=y2bar=z2bar=r2bar=0.0;
  for(ir=0;ir<NRMAX;ir++){
    r=(0.5+ir)*DELR;
    r2=r*r; 
    r3=r2*r;
    r4=r2*r2;
    if(!XSYM) xbar+=r3*A->GetElement(1,0,0,ir);
    if(!YSYM) ybar+=r3*A->GetElement(0,1,0,ir);
    if(!ZSYM) zbar+=r3*A->GetElement(0,0,1,ir);
    x2bar+=r4*A->GetElement(2,0,0,ir);
    y2bar+=r4*A->GetElement(0,2,0,ir);
    z2bar+=r4*A->GetElement(0,0,2,ir);
    r2bar+=r4*A->GetElement(0,0,0,ir);
  }
  xbar*=4.0*PI*DELR/3.0;
  ybar*=4.0*PI*DELR/3.0;
  zbar*=4.0*PI*DELR/3.0;
  printf("__________  EFFECTIVE GAUSSIAN PARAMETERS ____________\n");
  printf("Rinv=%g\n",sqrt(2.0*PI*DELR*r2bar/3.0));
  x2bar=4.0*PI*DELR*(2.0*x2bar/15.0+(r2bar/3.0))-xbar*xbar;
  y2bar=4.0*PI*DELR*(2.0*y2bar/15.0+(r2bar/3.0))-ybar*ybar;
  z2bar=4.0*PI*DELR*(2.0*z2bar/15.0+(r2bar/3.0))-zbar*zbar;
  printf("Gassian distribution with same offsets and 1-part. radii\n");
  printf("offset_xyz=(%g,%g,%g), R_xyz=(%g,%g,%g)\n",
	 xbar,ybar,zbar,
	 sqrt(fabs(0.5*x2bar)),sqrt(fabs(0.5*y2bar)),sqrt(fabs(0.5*z2bar)));
  printf("______________________________________________________\n");
}
void CSourceCalc::CalcEffGaussPars(CCHArray *A,double &Rx,double &Ry,
				   double &Rz,double &Xoff,double &Yoff,
				   double &Zoff){
  double xbar,ybar,zbar,x2bar,y2bar,z2bar,r2bar,r,r2,r3,r4,DELR;
  int ir,NRMAX;
  bool XSYM,YSYM,ZSYM;
  NRMAX=A->GetNRADIAL();
  DELR=A->GetRADSTEP();
  XSYM=A->GetXSYM();
  YSYM=A->GetYSYM();
  ZSYM=A->GetZSYM();
  const double PI=4.0*atan(1.0);
  xbar=ybar=zbar=x2bar=y2bar=z2bar=r2bar=0.0;
  for(ir=0;ir<NRMAX;ir++){
    r=(0.5+ir)*DELR;
    r2=r*r; 
    r3=r2*r;
    r4=r2*r2;
    if(!XSYM) xbar+=r3*A->GetElement(1,0,0,ir);
    if(!YSYM) ybar+=r3*A->GetElement(0,1,0,ir);
    if(!ZSYM) zbar+=r3*A->GetElement(0,0,1,ir);
    x2bar+=r4*A->GetElement(2,0,0,ir);
    y2bar+=r4*A->GetElement(0,2,0,ir);
    z2bar+=r4*A->GetElement(0,0,2,ir);
    r2bar+=r4*A->GetElement(0,0,0,ir);
  }
  xbar*=4.0*PI*DELR/3.0;
  ybar*=4.0*PI*DELR/3.0;
  zbar*=4.0*PI*DELR/3.0;
  printf("__________  EFFECTIVE GAUSSIAN PARAMETERS ____________\n");
  printf("Rinv=%g\n",sqrt(2.0*PI*DELR*r2bar/3.0));
  x2bar=4.0*PI*DELR*(2.0*x2bar/15.0+(r2bar/3.0))-xbar*xbar;
  y2bar=4.0*PI*DELR*(2.0*y2bar/15.0+(r2bar/3.0))-ybar*ybar;
  z2bar=4.0*PI*DELR*(2.0*z2bar/15.0+(r2bar/3.0))-zbar*zbar;
  printf("Gassian distribution with same offsets and 1-part. radii\n");
  printf("offset_xyz=(%g,%g,%g), R_xyz=(%g,%g,%g)\n",
	 xbar,ybar,zbar,
	 sqrt(fabs(0.5*x2bar)),sqrt(fabs(0.5*y2bar)),sqrt(fabs(0.5*z2bar)));
  printf("______________________________________________________\n");
  
  Xoff=xbar;
  Yoff=ybar;
  Zoff=zbar;
  Rx=sqrt(fabs(0.5*x2bar));
  Ry=sqrt(fabs(0.5*y2bar));
  Rz=sqrt(fabs(0.5*z2bar));

}

void CSourceCalc::NormCheck(C3DArray *threed){
  int nsx,nsy,nsz,isx,isy,isz,ix,iy,iz;
  int nxmax=threed->GetNXMAX();
  int nymax=threed->GetNYMAX();
  int nzmax=threed->GetNZMAX();
  double prefactor=threed->GetDELX()*threed->GetDELY()*threed->GetDELZ();
  double norm=0.0;
  nsx=nsy=nsz=2;
  if(threed->GetXSYM()) nsx=1;
  if(threed->GetYSYM()) nsy=1;
  if(threed->GetZSYM()) nsz=1;
  if(nsx==1) prefactor*=2.0;
  if(nsy==1) prefactor*=2.0;
  if(nsz==1) prefactor*=2.0;
  for(ix=0;ix<nxmax;ix++){
    for(iy=0;iy<nymax;iy++){
      for(iz=0;iz<nzmax;iz++){
	for(isz=0;isz<nsz;isz++){
	  for(isy=0;isy<nsy;isy++){
	    for(isx=0;isx<nsx;isx++){
	      norm+=threed->GetElement(isx,ix,isy,iy,isz,iz)*prefactor;
	    }
	  }
	}
      }
    }
  }
  printf("Norm Check of 3DArray = %g\n",norm);
}


#include "sourcecalc_blast.cc"
#include "sourcecalc_gauss.cc"
#include "sourcecalc_OSCAR.cc"
#include "sourcecalc_GX1d.cc"


#endif
