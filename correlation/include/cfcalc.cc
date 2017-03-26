#ifndef __INCLUDE_CFCalc_CC__
#define __INCLUDE_CFCalc_CC__

#include "cfcalc.h"
using namespace std;

double CFCalc::GetChiSquared(C3DArray *CFexp,C3DArray *Error,
			       C3DArray *CFtheory){
  if(!ArrayCalc::CompareArrayParameters(CFexp,Error)) exit(1);
  if(!ArrayCalc::CompareArrayParameters(CFexp,CFtheory)) exit(1);
  double chisquared=0.0,numer,denom;
  int ix,iy,iz,isx,isy,isz,nsx,nsy,nsz,nxmax,nymax,nzmax;
  int ndof=0;
  nsx=nsy=nsz=2;
  if(CFexp->GetXSYM()) nsx=1;
  if(CFexp->GetYSYM()) nsy=1;
  if(CFexp->GetZSYM()) nsz=1;
  nxmax=CFexp->GetNXMAX();
  nymax=CFexp->GetNYMAX();
  nzmax=CFexp->GetNZMAX();
  for(isx=0;isx<nsx;isx++){
    for(ix=0;ix<nxmax;ix++){
      for(isy=0;isy<nsy;isy++){
	for(iy=0;iy<nymax;iy++){
	  for(isz=0;isz<nsz;isz++){
	    for(iz=0;iz<nzmax;iz++){
	      numer=pow(CFexp->GetElement(isx,ix,isy,iy,isz,iz)
			-CFtheory->GetElement(isx,ix,isy,iy,isz,iz),2);
	      denom=Error->GetElement(isx,ix,isy,iy,isz,iz);
	      if(denom<1.0E-8) printf("CFCalc::GetChiSquared : \nSuspiciously small or negative error = %g, isx=%d,ix=%d, isy=%d,iy=%d, isz=%d,iz=%d\n",
				      denom,isx,ix,isy,iy,isz,iz);
	      denom=denom*denom;
	      chisquared+=numer/denom;
	      ndof+=1;
	    }
	  }
	}
      }
    }
  }
  printf("chi^2/Ndof=%g\n",chisquared/double(ndof));
  return chisquared;
}

double CFCalc::GetChiSquared(int lx,int ly,int lz,CCHArray *CFexp,
			     CCHArray *Error,CCHArray *CFtheory){
  int ir,nrmax;
  double chisquared=0.0,numer,denom;
  nrmax=CFexp->GetNRADIAL();
  for(ir=0;ir<nrmax;ir++){
    denom=Error->GetElement(lx,ly,lz,ir);
    if(denom<1.0E-8) printf("CFCalc::GetChiSquared : \nSuspiciously small or negative error = %g, l=(%,%d,%d)d\n",denom,lx,ly,lz);
    numer=CFexp->GetElement(lx,ly,lz,ir)
      -CFtheory->GetElement(lx,ly,lz,ir);
    chisquared+=numer*numer/(denom*denom);
  }
  printf("chi^2/Ndof=%g\n",chisquared/double(nrmax));
  return chisquared;
}

#endif
