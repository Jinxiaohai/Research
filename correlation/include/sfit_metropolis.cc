#ifndef __INCLUDE_SFIT_METROPOLIS_CC__
#define __INCLUDE_SFIT_METROPOLIS_CC__

#include "sfit_metropolis.h"

void CParInfo::Set(string nameset,double xset,double errorset,
		double xminset,double xmaxset){
  strcpy(name,nameset.c_str());
  x=xset;
  error=errorset;
  xmin=xminset;
  xmax=xmaxset;
}

void CCF2S_Metropolis::PrintPars(){
  int ipar;
  
  printf("ipar        name     value       error       min         max     fixed\n");
  for(ipar=0;ipar<npars;ipar++){
    printf("%2d : %12s %11.4e %11.4e %11.4e %11.4e  %d\n",
	   ipar,par[ipar].name,par[ipar].x,par[ipar].error,
	   par[ipar].xmin,par[ipar].xmax,par[ipar].fixed);
  }
}

void CCF2S_Metropolis::Init(){
  int i;
  chisquared_best=1.0E20;
  nfreepars=npars;
  randy=new CRandom(-1234);
  par=new CParInfo[npars];
  x_best=new double[npars];
  matrix=new CGSLMatrix_Real(npars);
  sigma=new double *[npars];
  Usigma=new double *[npars];
  xran=new double[npars];
  Ueigen=new double[npars];
  oldchisquared=1.0E10;
  chisquared=oldchisquared;
  for(i=0;i<npars;i++){
    par[i].x=par[i].oldx=par[i].xmin=par[i].xmax=0.0;
    par[i].xmin=-1.0E10; par[i].xmax=1.0E10;
    sigma[i]=new double[npars];
    Usigma[i]=new double[npars];
    par[i].fixed=0;
  }
  Reset();
  printf("Metropolis Arrays Initialized\n");
}

void CCF2S_Metropolis::FixPar(int ipar){
  nfreepars-=1;
  par[ipar].fixed=1;
}

void CCF2S_Metropolis::FreePar(int ipar){
  nfreepars+=1;
  par[ipar].fixed=0;
}

void CCF2S_Metropolis::Minimize(int maxcalls){
  Reset();
  int icall,ipar,j;
  double xarg,stepscale,step;
  stepscale=1.0/sqrt(double(nfreepars));
  for(icall=0;icall<maxcalls;icall++){
    oldchisquared=chisquared;
    for(ipar=0;ipar<npars;ipar++) par[ipar].oldx=par[ipar].x;
    for(ipar=0;ipar<npars;ipar++){
      if(par[ipar].fixed==0){
	do{
	  par[ipar].x=par[ipar].oldx;
	  for(j=0;j<npars;j++) xran[j]=randy->gauss();
	  for(j=0;j<npars;j++){
	    step=Usigma[ipar][j]*xran[j]*stepscale;
	    //printf("U(%d,%d)=%g, step=%g \n",ipar,j,Usigma[ipar][j],step);
	    par[ipar].x+=step;
	  }
	} while(par[ipar].x<par[ipar].xmin || par[ipar].x>par[ipar].xmax);
      }
    }
    CalcChiSquare();
    if(chisquared>oldchisquared){
      xarg=0.5*(chisquared-oldchisquared);
      if(xarg>20 || randy->ran()>exp(-xarg)){
	chisquared=oldchisquared;
	for(ipar=0;ipar<npars;ipar++) par[ipar].x=par[ipar].oldx;
      }
      else Nsuccess+=1;
    }
    else Nsuccess+=1;
    
    for(ipar=0;ipar<npars;ipar++){
      par[ipar].xbar+=par[ipar].x;
      for(j=0;j<npars;j++) sigma[ipar][j]+=par[ipar].x*par[j].x;
    }
    if(chisquared<chisquared_best){
      chisquared_best=chisquared;
      for(ipar=0;ipar<npars;ipar++) x_best[ipar]=par[ipar].x;
    }
    printf("nsuccess=%d\n",Nsuccess);

  }
  printf("_____ Best chi^2=%g _____\nBest x[] = ",chisquared_best);
  for(ipar=0;ipar<npars;ipar++) printf("%g ",x_best[ipar]);
  printf("\n");
  UpdateSigma();
}
  
void CCF2S_Metropolis::Reset(){
  int i,j;
  Ncalls=0;
  Nsuccess=0;
  for(i=0;i<npars;i++){
    par[i].xbar=0.0;
    for(j=0;j<npars;j++) sigma[i][j]=0.0;
  }
}

void CCF2S_Metropolis::UpdateStepSize(){
  int i,j;
  matrix->EigenFind(sigma,Usigma,Ueigen);
  for(i=0;i<npars;i++){
    if(Ueigen[i]<-1.0E-10){
      printf("FATAL: In UpdateStepSize, negative eigenvalue, =%g\n",Ueigen[i]);
      exit(1);
    }
  }
  for(i=0;i<npars;i++){
    par[i].error=0.0;
    for(j=0;j<npars;j++){
      Usigma[i][j]=Usigma[i][j]*sqrt(fabs(Ueigen[j]));
      par[i].error+=Usigma[i][j]*Usigma[i][j];
    }
    par[i].error=sqrt(par[i].error);
  }
}

void CCF2S_Metropolis::UpdateSigma(){
  int i,j;
  for(i=0;i<npars;i++){
    par[i].xbar=par[i].xbar/double(Ncalls);
      for(j=0;j<npars;j++) sigma[i][j]=sigma[i][j]/double(Ncalls);
  }
  for(i=0;i<npars;i++)
    for(j=0;j<npars;j++) sigma[i][j]=(sigma[i][j]-par[i].xbar*par[j].xbar)
			   *double(Nsuccess)/double(Nsuccess-1);
}


void CCF2S_Metropolis::PrintSigma(){
  int ipara,iparb;
  printf("________ < delX_i delX_j > ________\n");
  for(ipara=0;ipara<npars;ipara++){
    for(iparb=0;iparb<npars;iparb++) printf("%9.2e ",sigma[ipara][iparb]);
    printf("\n");
  }
}

void CCF2S_Metropolis::PrintUsigma(){
  int ipara,iparb;
  printf("________ < USigma_ij > ________\n");
  for(ipara=0;ipara<npars;ipara++){
    for(iparb=0;iparb<npars;iparb++) printf("%9.2e ",Usigma[ipara][iparb]);
    printf("\n");
  }
}

CCF2S_Metropolis::CCF2S_Metropolis(){
}


void CCF2S_Metropolis::CalcChiSquare(){
  int ipar;

  Ncalls+=1;
  printf("In CalcChiSquare, x[] = ");
  for(ipar=0;ipar<npars;ipar++) printf("%g,",par[ipar].x);
  printf("\n");
  for(ipar=0;ipar<npars;ipar++){
    parameter::set(sourcecalc->spars,par[ipar].name,par[ipar].x);
  }
  if(ndim==1){
    sourcecalc->CalcS(lx,ly,lz,source);
    S2CF::s2c(source,kernel,ctheory);
    chisquared=CFCalc::GetChiSquared(lx,ly,lz,cexp,cerror,ctheory);
  }
  else if(ndim==3){
    sourcecalc->CalcS(source);
    S2CF::s2c(source,kernel,ctheory);
    ArrayCalc::Calc3DArrayFromAExpArray(ctheory,ctheory3D);
    chisquared=CFCalc::GetChiSquared(cexp3D,cerror3D,ctheory3D);
  }
  printf("chi^2=%g   (Ncalls=%d)\n",chisquared,Ncalls);
}

#include "sfit_metropolis_3dgaussian.cc"
//#include "sfit_metropolis_GX1d.cc"
#include "sfit_metropolis_blast.cc"


#endif

