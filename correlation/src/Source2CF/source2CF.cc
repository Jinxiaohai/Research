#ifndef __INCLUDE_S2C_CC
#define __INCLUDE_S2C_CC

#include "source2CF.h"

using namespace std;

void S2CF::s2c(int lx,int ly,int lz,CCHArray *S,CKernel *kernel,CCHArray *CF){
  int lmax=kernel->GetLMAX();
  if(lmax>CF->GetLMAX() || lmax>S->GetLMAX()){
    printf("FATAL: Kernel LMAX=%d larger than either S LMAX=%d or CF LMAX=%d\n",
	   lmax,S->GetLMAX(),CF->GetLMAX());
    exit(1);
  }
  if( (CF->GetLMAX()!=lmax) ){
    printf("WARNING: Array parameters for kernel calculations don't match\n");
    printf("___ CORRELATION FUNCTION PARAMETERS ___\n");
    CF->PrintPars();
    printf("_____ SOURCE FUNCTION PARAMETERS _____\n");
    S->PrintPars();
    printf("For kernel, LMAX=%d\n",kernel->GetLMAX());
  }
  if(kernel->GetIDENTICAL() && (!CF->GetXSYM() || !CF->GetYSYM() || !CF->GetZSYM())){
    printf("FATAL: kernel has no odd L components, but CF wants them\n");
    printf("Make sure CF array has XSYM=YSYM=ZSYM='true'\n");
    exit(1);
  }
  int iq,ir,nqmax,nrmax,L;
  double r,delr,q,delq,norm;
  const double PI=4.0*atan(1.0);
  bool match=0;
  delr=S->GetRADSTEP();
  nrmax=S->GetNRADIAL();
  delq=CF->GetRADSTEP();
  nqmax=CF->GetNRADIAL();
  if(fabs(delr-kernel->GetDELR())<1.0E-5 
     && fabs(delq-kernel->GetDELQ())<1.0E-5
     && nrmax==kernel->GetNRMAX()
     && nqmax==kernel->GetNQMAX()) match=1;
  CF->ZeroArray(lx,ly,lz);
  L=lx+ly+lz;
  for(iq=0;iq<nqmax;iq++){
    q=(0.5+iq)*delq;
    for(ir=0;ir<nrmax;ir++){
      r=(0.5+ir)*delr;
      norm=4.0*PI*r*r*delr;
      if(match) CF->IncrementElement(lx,ly,lz,iq,kernel->GetValue(L,iq,ir)
				     *norm*S->GetElement(lx,ly,lz,ir));
      else CF->IncrementElement(lx,ly,lz,iq,kernel->GetValue(L,q,r)
				*norm*S->GetElement(lx,ly,lz,r));
    }
  }
}

void S2CF::s2c(CCHArray *S,CKernel *kernel,CCHArray *CF){
  int lmax=kernel->GetLMAX();
  if(lmax>CF->GetLMAX() || lmax>S->GetLMAX()){
    printf("FATAL: Kernel LMAX=%d larger than either S LMAX=%d or CF LMAX=%d\n",
	   lmax,S->GetLMAX(),CF->GetLMAX());
    exit(1);
  }
  if( (CF->GetLMAX()!=lmax) ){
    printf("WARNING: Array parameters for kernel calculations don't match\n");
    printf("___ CORRELATION FUNCTION PARAMETERS ___\n");
    CF->PrintPars();
    printf("_____ SOURCE FUNCTION PARAMETERS _____\n");
    S->PrintPars();
    printf("For kernel, LMAX=%d\n",kernel->GetLMAX());
  }
  if(kernel->GetIDENTICAL() && (!CF->GetXSYM() || !CF->GetYSYM() || !CF->GetZSYM())){
    printf("FATAL: kernel has no odd L components, but CF wants them\n");
    printf("Make sure CF array has XSYM=YSYM=ZSYM='true'\n");
    exit(1);
  }
  int iq,ir,nqmax,nrmax,L,lx,ly,lz,dlx,dly,dlz;
  double r,delr,q,delq,norm;
  const double PI=4.0*atan(1.0);
  bool match=0;
  delr=S->GetRADSTEP();
  nrmax=S->GetNRADIAL();
  delq=CF->GetRADSTEP();
  nqmax=CF->GetNRADIAL();
  dlx=dly=dlz=1;
  if(CF->GetXSYM()) dlx=2;
  if(CF->GetYSYM()) dly=2;
  if(CF->GetZSYM()) dlz=2;
  if(fabs(delr-kernel->GetDELR())<1.0E-5 
     && fabs(delq-kernel->GetDELQ())<1.0E-5
     && nrmax==kernel->GetNRMAX()
     && nqmax==kernel->GetNQMAX()) match=1;
  CF->ZeroArray();
  for(iq=0;iq<nqmax;iq++){
    q=(0.5+iq)*delq;
    for(lx=0;lx<2;lx+=dlx){
      for(ly=0;ly<=lmax-lx;ly+=dly){
	for(lz=0;lz<=lmax-lx-ly;lz+=dlz){
	  L=lx+ly+lz;
	  for(ir=0;ir<nrmax;ir++){
	    r=(0.5+ir)*delr;
	    norm=4.0*PI*r*r*delr;
	    if(match) CF->IncrementElement(lx,ly,lz,iq,kernel->GetValue(L,iq,ir)
					   *norm*S->GetElement(lx,ly,lz,ir));
	    else CF->IncrementElement(lx,ly,lz,iq,kernel->GetValue(L,q,r)
				      *norm*S->GetElement(lx,ly,lz,r));
	  }
	}
      }
    }
  }
}

void S2CF::s2c(C3DArray *S,CKernelWF *kernel,C3DArray *CF){
  int ix,iy,iz,isx,isy,isz,jx,jy,jz,jsx,jsy,jsz;
  int nsx,nsy,nsz;
  int nxmax=S->GetNXMAX();
  int nymax=S->GetNYMAX();
  int nzmax=S->GetNZMAX();
  double delx=S->GetDELX();
  double dely=S->GetDELY();
  double delz=S->GetDELZ();
  int nqxmax=CF->GetNXMAX();
  int nqymax=CF->GetNYMAX();
  int nqzmax=CF->GetNZMAX();
  double delqx=CF->GetDELX();
  double delqy=CF->GetDELY();
  double delqz=CF->GetDELZ();
  double x,y,z,qx,qy,qz,norm,q,r,ctheta,wf2,svalue;

  norm=delx*dely*delz;
  nsx=nsy=nsz=2;
  if(S->GetXSYM()){
    nsx=1;
    norm*=2.0;
  }
  if(S->GetYSYM()){
    nsy=1;
    norm*=2.0;
  }
  if(S->GetZSYM()){
    nsz=1;
    norm*=2.0;
  }
  CF->ZeroArray();

  printf("beginning big loop in s2c\n");
  printf("nsxyz=(%d,%d,%d), nxyzmax=(%d,%d,%d), nqxyzmax=(%d,%d,%d)\n",
	 nsx,nsy,nsz,nxmax,nymax,nzmax,nqxmax,nqymax,nqzmax);
  for(isx=0;isx<nsx;isx++){
    for(ix=0;ix<nxmax;ix++){
      x=(0.5+ix)*delx;
      if(isx==1) x=-x;
      for(isy=0;isy<nsy;isy++){
	for(iy=0;iy<nymax;iy++){
	  y=(0.5+iy)*dely;
	  if(isy==1) y=-y;
	  for(isz=0;isz<nsz;isz++){
	    for(iz=0;iz<nzmax;iz++){
	      z=(0.5+iz)*delz;
	      if(isz==1) z=-z;
	      r=sqrt(x*x+y*y+z*z);
	      svalue=S->GetElement(isx,ix,isy,iy,isz,iz);
	      for(jsx=0;jsx<nsx;jsx++){
		for(jx=0;jx<nqxmax;jx++){
		  qx=(0.5+jx)*delqx;
		  if(jsx==1) qx=-qx;
		  for(jsy=0;jsy<nsy;jsy++){
		    for(jy=0;jy<nqymax;jy++){
		      qy=(0.5+jy)*delqy;
		      if(jsy==1) qy=-qy;
		      for(jsz=0;jsz<nsz;jsz++){
			for(jz=0;jz<nqzmax;jz++){
			  qz=(0.5+jz)*delqz;
			  if(jsz==1) qz=-qz;
			  q=sqrt(qx*qx+qy*qy+qz*qz);
			  ctheta=(qx*x+qy*y+qz*z)/(q*r);
			  wf2=kernel->GetPsiSquared(q,r,ctheta);
			  //printf("ixyz=(%d,%d,%d), jxyz=(%d,%d,%d)\n",
			  //	 ix,iy,iz,jx,jy,jz);
			  //printf("q=%g,r=%g,ctheta=%g,wf2=%g\n",q,r,ctheta,wf2);
			  CF->IncrementElement(jsx,jx,jsy,jy,jsz,jz,norm*wf2*svalue);
			  
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  printf("ending big loop in s2c\n");

}


#endif
