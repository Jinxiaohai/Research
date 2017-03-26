#ifndef __INCLUDE_SFIT_CC__
#define __INCLUDE_SFIT_CC__

#include "sfit.h"

void CParInfo::Set(string nameset,double xset,double errorset,
		double xminset,double xmaxset){
  strcpy(name,nameset.c_str());
  currentx=xset;
  error=errorset;
  xmin=xminset;
  xmax=xmaxset;
}
CParInfo::CParInfo(){
  name=new char[20];
}
CParInfo::~CParInfo(){
  delete [] name;
}

void CCF2SFit::SetPar(string parstring,double xset){
  int i;
  char parname[20];
  strcpy(parname,parstring.c_str());
  i=0;
  while(Misc::comparestrings(par[i]->name,parname)==0 && i<npars){
    i+=1;
  }
  if(i<npars) par[i]->currentx=xset;
  else printf("Can not set %s, parameter with that name does not exist\n",
	      parname);
}

void CCF2SFit::UseBestPars(){
  int i;
  currentchisquared=bestchisquared;
  for(i=0;i<npars;i++) par[i]->currentx=par[i]->bestx;
}

void CCF2SFit::SetL(int lxset,int lyset,int lzset){
  lx=lxset;
  ly=lyset;
  lz=lzset;
}

void CCF2SFit::SwitchPars(int i,int j){
  CParInfo *partemp;
  partemp=par[j];
  par[j]=par[i];
  par[i]=partemp;

}

void CCF2SFit::FreePar(string parstring){
  int i;
  char parname[20];
  strcpy(parname,parstring.c_str());
  i=nfreepars;
  while(Misc::comparestrings(par[i]->name,parname)==0 && i<npars){
    i+=1;
  }
  if(i==npars) printf("Par %s already free or does not exist\n",
		      parstring.c_str());
  else{
    if(i!=nfreepars) SwitchPars(i,nfreepars);
    par[nfreepars]->fixed=0;
    printf("Freeing par[%d], name=%s\n",nfreepars,par[nfreepars]->name);
    nfreepars+=1;
  }
}

void CCF2SFit::FixPar(string parstring){
  int i;
  char parname[20];
  strcpy(parname,parstring.c_str());
  i=0;
  while(Misc::comparestrings(par[i]->name,parname)==0 && i<nfreepars){
    i+=1;
  }
  if(i==nfreepars) printf("Par %s already fixed or does not exist\n",
		      parstring.c_str());
  else{
    if(i!=nfreepars-1) SwitchPars(i,nfreepars-1);
    nfreepars-=1;
    par[nfreepars]->fixed=1;
    printf("Fixing par[%d], name=%s\n",nfreepars,par[nfreepars]->name);
  }
}

void CCF2SFit::PrintPars(){
  int i;
  printf("ipar        name     value       error       min         max     fixed\n");
  for(i=0;i<npars;i++){
    printf("%2d : %12s %11.4e %11.4e %11.4e %11.4e  %d\n",
	   i,par[i]->name,par[i]->currentx,par[i]->error,
	   par[i]->xmin,par[i]->xmax,par[i]->fixed);
  }
}

void CCF2SFit::Init(){
  int i,j;
  currentchisquared=1.0E20;
  bestchisquared=1.0E20;
  nfreepars=npars;
  ncalls=0;
  randy=new CRandom(-1234);
  par=new CParInfo *[npars];
  ErrorMatrix=new double *[npars];
  StepMatrix=new double *[npars];
  for(i=0;i<npars;i++){
    par[i]=new CParInfo();
    par[i]->bestx=0.0;
    par[i]->xmin=-1.0E10; par[i]->xmax=1.0E10;
    ErrorMatrix[i]=new double[npars];
    StepMatrix[i]=new double[npars];
    par[i]->fixed=0;
    for(j=0;j<npars;j++){
      ErrorMatrix[i][j]=0.0;
      StepMatrix[i][j]=0.0;
    }
  }
  printf("SFit Arrays Initialized\n");
}

void CCF2SFit::SteepestDescent(int maxtries){
  int i,j,itry,nfailure,nsuccess;
  double chisquared,newchisquared,chi2a,chi2b,qstep,maxqstep=3.0;
  double d2chi2dq2,dchi2dq_mag,dq;
  double *delx,*xstep,*delq,*x,*xa,*xb,*qxratio,*xnew;
  double *dchi2dq;
  delx=new double[npars];
  x=new double[npars];
  xa=new double[npars];
  xb=new double[npars];
  xnew=new double[npars];
  delq=new double[npars];
  qxratio=new double[npars];
  dchi2dq=new double[npars];
  for(i=0;i<npars;i++){
    qxratio[i]=par[i]->error;
    if(par[i]->fixed) qxratio[i]=0.0;
    x[i]=par[i]->currentx;
    dchi2dq[i]=delx[i]=0.0;
  }
  dq=3.0;
  chisquared=GetChiSquared(x);

  nsuccess=0;
  itry=0;
  do{
    itry+=1;
    dchi2dq_mag=0.0;
    for(i=0;i<npars;i++){
      if(par[i]->fixed==0){
	delx[i]=qxratio[i]*dq;
	for(j=0;j<npars;j++){
	  xa[j]=x[j];
	  xb[j]=x[j];
	}
	xa[i]=x[i]-0.5*delx[i];
	xb[i]=x[i]+0.5*delx[i];
	chi2a=GetChiSquared(xa);
	chi2b=GetChiSquared(xb);
	dchi2dq[i]=(chi2b-chi2a)*qxratio[i]/delx[i];
	dchi2dq_mag+=dchi2dq[i]*dchi2dq[i];
      }
    }
    dchi2dq_mag=sqrt(dchi2dq_mag);

    for(i=0;i<npars;i++){
      if(par[i]->fixed==0){
	delq[i]=-dq*dchi2dq[i]/dchi2dq_mag;
	delx[i]=delq[i]*qxratio[i];
	xa[i]=x[i]-0.5*delx[i];
	xb[i]=x[i]+0.5*delx[i];
	printf("xa[%d]=%g, xb[%d]=%g\n",i,xa[i],i,xb[i]);
      }
    }
    chi2a=GetChiSquared(xa);
    chi2b=GetChiSquared(xb);
    d2chi2dq2=4.0*(chi2b+chi2a-2.0*chisquared)/(dq*dq);
    qstep=dchi2dq_mag/d2chi2dq2;
    if(qstep<0.0){
      qstep=dq;
      printf("Warning: curvature in q has wrong sign\n");
    }
    if(qstep>10.0*dq) qstep=10.0*dq;
  TRYNEWQSTEP:
    nfailure=0;
    for(i=0;i<npars;i++){
      if(par[i]->fixed==0){
	delq[i]=-qstep*dchi2dq[i]/dchi2dq_mag;
	xnew[i]=x[i]+delq[i]*qxratio[i];
      }
    }
    newchisquared=GetChiSquared(xnew);
    if(newchisquared>chisquared){
      printf("new chi^2 not so good\n");
      nfailure+=1;
      if(nfailure>5){
	printf("STEEPEST DESCENT FAILURE\n");
	exit(1);
      }
      qstep=0.5*qstep;
      goto TRYNEWQSTEP;
    }
    if(qstep<2.0*dq){
      nsuccess+=1;
      dq=dq/3.0;
      chisquared=newchisquared;
    }

    for(i=0;i<npars;i++){
      if(par[i]->fixed==0) x[i]=xnew[i];
    }
  }while(itry<maxtries && nsuccess<2);
  
  chisquared=GetChiSquared(x);
  currentchisquared=chisquared;
  for(i=0;i<npars;i++){
    if(par[i]->fixed==0) par[i]->currentx=x[i];
  }

  delete [] xa;
  delete [] xb;
  delete [] x;
  delete [] xnew;
  delete [] delx;
  delete [] delq;
  delete [] dchi2dq;
  delete [] qxratio;
}

void CCF2SFit::Metropolis(int maxcalls){
  int icall,i,j,Nsuccess=0;
  bool success;
  double stepscale,step,chisquared;
  double *x,*xran,*xbar;
  x=new double[npars];
  xran=new double[npars];
  xbar=new double[npars];
  for(i=0;i<npars;i++){
    x[i]=par[i]->currentx;
    xbar[i]=0.0;
    for(j=0;j<npars;j++) ErrorMatrix[i][j]=0.0;
  }
  stepscale=1.0/sqrt(double(nfreepars));
  for(icall=0;icall<maxcalls;icall++){

  GETNEWXRAN:
    for(j=0;j<npars;j++) xran[j]=randy->gauss();
    for(i=0;i<npars;i++){
      if(par[i]->fixed==0){
	x[i]=par[i]->currentx;
	for(j=0;j<npars;j++){
	  step=StepMatrix[i][j]*xran[j]*stepscale;
	  //printf("U(%d,%d)=%g, step=%g \n",i,j,StepMatrix[i][j],step);
	  x[i]+=step;
	}
      }
    }
    for(i=0;i<npars;i++)
      if(x[i]<par[i]->xmin || x[i]>par[i]->xmax) goto GETNEWXRAN;

    chisquared=GetChiSquared(x);
    success=0;
    if(chisquared<currentchisquared){
      success=1;
    }
    else if(randy->ran()<exp(-0.5*(chisquared-currentchisquared))){
      success=1;
    }

    if(success==1){
      Nsuccess+=1;
      currentchisquared=chisquared;
      for(i=0;i<npars;i++) par[i]->currentx=x[i];
      printf("SUCCESS, Nsuccess=%d, Ncalls=%d\n",Nsuccess,icall+1);
    }

    for(i=0;i<npars;i++){
      xbar[i]+=x[i];
      for(j=0;j<npars;j++) ErrorMatrix[i][j]+=x[i]*x[j];
    }
  }

  for(i=0;i<npars;i++){
    par[i]->xbar=xbar[i]/double(maxcalls);
    for(j=0;j<npars;j++) ErrorMatrix[i][j]=ErrorMatrix[i][j]/double(maxcalls);
  }
  for(i=0;i<npars;i++)
    for(j=0;j<npars;j++) ErrorMatrix[i][j]=(ErrorMatrix[i][j]-par[i]->xbar*par[j]->xbar)
			   *double(Nsuccess)/double(Nsuccess-1);
  
  printf("Best chi^2=%g, Best x[] = ",bestchisquared);
  for(i=0;i<npars;i++) printf("%g ",par[i]->bestx);
  printf("\n");

  delete [] x;
  delete [] xran;
  delete [] xbar;
}
  
void CCF2SFit::UpdateMetroStepSize(){
  int i,j;
  double *SMeigenval;
  CGSLMatrix_Real *matrixcalc;
  matrixcalc=new CGSLMatrix_Real(npars);
  SMeigenval=new double[npars];
  matrixcalc->EigenFind(ErrorMatrix,StepMatrix,SMeigenval);
  for(i=0;i<npars;i++){
    if(SMeigenval[i]<-1.0E-10){
      printf("FATAL: In UpdateStepSize, negative eigenvalue, =%g\n",
	     SMeigenval[i]);
      exit(1);
    }
  }
  for(i=0;i<npars;i++){
    par[i]->error=0.0;
    for(j=0;j<npars;j++){
      StepMatrix[i][j]=StepMatrix[i][j]*sqrt(fabs(SMeigenval[j]));
      par[i]->error+=StepMatrix[i][j]*StepMatrix[i][j];
    }
    par[i]->error=sqrt(par[i]->error);
  }
  delete [] SMeigenval;
  delete(matrixcalc);
}


void CCF2SFit::PrintErrorMatrix(){
  int ia,ib;
  printf("________ < delX_i delX_j > ________\n");
  for(ia=0;ia<npars;ia++){
    for(ib=0;ib<npars;ib++) printf("%9.2e ",ErrorMatrix[ia][ib]);
    printf("\n");
  }
  printf("___________________________________\n");
}

void CCF2SFit::PrintStepMatrix(){
  int ia,ib;
  printf("________ < StepMatrix_ij > ________\n");
  for(ia=0;ia<npars;ia++){
    for(ib=0;ib<npars;ib++) printf("%9.2e ",StepMatrix[ia][ib]);
    printf("\n");
  }
  printf("___________________________________\n");
}

CCF2SFit::CCF2SFit(){
}
CCF2SFit::~CCF2SFit(){
  int i;
  delete [] par;
  for(i=0;i<npars;i++){
    delete [] ErrorMatrix[i];
    delete [] StepMatrix[i];
  }
  delete [] ErrorMatrix;
  delete [] StepMatrix;
}


double CCF2SFit::GetChiSquared(double *xx){
  double chisquared;
  int i;
  double *x;
  ncalls+=1;
  x=new double[npars];
  for(i=0;i<nfreepars;i++) x[i]=xx[i];
  for(i=nfreepars;i<npars;i++) x[i]=par[i]->currentx;

  printf("In GetChiSquare, x[] = ");
  for(i=0;i<npars;i++) printf("%g,",x[i]);
  printf("\n");
  for(i=0;i<npars;i++){
    parameter::set(sourcecalc->spars,par[i]->name,x[i]);
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
  if(chisquared<bestchisquared){
    bestchisquared=chisquared;
    for(i=0;i<npars;i++) par[i]->bestx=x[i];
  }
  return chisquared;
}

void CCF2SFit::Newton(int maxtries){
  bool success;
  int i,j,k,itry,nrescale,nsuccess;
  double **curvature,*slope;
  double *x,*xnew,*xa,*xb,*xc,*xd,*delx,*dx;
  double chi2,chi2a,chi2b,chi2c,chi2d,newchi2,scheck;
  CGSLMatrix_Real *matrixcalc;
  matrixcalc=new CGSLMatrix_Real(nfreepars);
  slope=new double[nfreepars];
  x=new double[nfreepars];
  xnew=new double[nfreepars];
  xa=new double[nfreepars];
  xb=new double[nfreepars];
  xc=new double[nfreepars];
  xd=new double[nfreepars];
  delx=new double[nfreepars];
  dx=new double[nfreepars];
  curvature=new double*[nfreepars];
  for(i=0;i<nfreepars;i++) curvature[i]=new double[nfreepars];

  for(i=0;i<nfreepars;i++){
    dx[i]=par[i]->error;
    x[i]=par[i]->currentx;
  }

  nsuccess=itry=0;
  do{
    itry+=1;
    chi2=GetChiSquared(x);
    for(i=1;i<nfreepars;i++){
      for(j=0;j<i;j++){
	for(k=0;k<nfreepars;k++) xa[k]=xb[k]=xc[k]=xd[k]=x[k];
	xa[i]+=dx[i]; xa[j]+=dx[j];
	xb[i]+=dx[i]; xb[j]-=dx[j];
	xc[i]-=dx[i]; xc[j]+=dx[j];
	xd[i]-=dx[i]; xd[j]-=dx[j];
	chi2a=GetChiSquared(xa);
	chi2b=GetChiSquared(xb);
	chi2c=GetChiSquared(xc);
	chi2d=GetChiSquared(xd);
	curvature[i][j]=(chi2a+chi2d-chi2b-chi2c)/(4.0*dx[i]*dx[j]);
	curvature[j][i]=curvature[i][j];
	printf("___ itry=%d, Curvature[%d][%d]=%g ___\n",
	       itry,i,j,curvature[i][j]);
      }
    }
    for(i=0;i<nfreepars;i++){
      for(k=0;k<nfreepars;k++) xa[k]=xb[k]=x[k];
      xa[i]+=dx[i];
      xb[i]-=dx[i];
      chi2a=GetChiSquared(xa);
      chi2b=GetChiSquared(xb);
      curvature[i][i]=(chi2a-2*chi2+chi2b)/(dx[i]*dx[i]);
      printf("___ itry=%d, Curvature[%d][%d]=%g ___\n",itry,i,i,curvature[i][i]);
      slope[i]=(chi2a-chi2b)/(2.0*dx[i]);
    }
    matrixcalc->SolveLinearEqs(slope,curvature,delx);
    nrescale=0;
  TRY_RESCALED:
    printf("delx[] = ");
    for(i=0;i<nfreepars;i++){
      printf("%g ",delx[i]);
      xnew[i]=x[i]-delx[i];
    }
    printf("\n");

    for(i=0;i<nfreepars;i++){
      if(par[i]->fixed==0){
	if(xnew[i]<par[i]->xmin || xnew[i]>par[i]->xmax){
	  for(j=0;j<nfreepars;j++) delx[j]=0.5*delx[j];
	  printf("Stepped outside min/max, will rescale delx[]\n");
	  if(nrescale>4) exit(1);
	  nrescale+=1;
	  goto TRY_RESCALED;
	}
      }
    }
    
    newchi2=GetChiSquared(xnew);
    if(newchi2>chi2){
      for(j=0;j<nfreepars;j++){
	delx[j]=0.5*delx[j];
      }
      printf("new chi^2 bigger than previous, will rescale delx[]\n");
      goto TRY_RESCALED;
    }
    
    success=1;
    scheck=0.0;
    for(i=0;i<nfreepars;i++){
      for(j=0;j<nfreepars;j++){
	scheck+=curvature[i][j]*delx[i]*delx[j];
      }
    }
    if(scheck<1.0) success=1;

    if(success){
      printf("SUCCESS!!!!!!!!!!!\n");
      nsuccess+=1;
      CalcErrorMatrixFromCurvature(curvature);
      for(i=0;i<nfreepars;i++) dx[i]=par[i]->error;
    }

    chi2=newchi2;
    for(i=0;i<nfreepars;i++) x[i]=xnew[i];
  } while(itry<maxtries && nsuccess<2);

  currentchisquared=chi2;
  for(i=0;i<nfreepars;i++) par[i]->currentx=x[i];

  delete [] dx;
  delete [] x;
  delete [] xnew;
  delete [] delx;
  delete [] xa;
  delete [] xb;
  delete [] xc;
  delete [] xd;
  delete [] slope;
  for(i=0;i<nfreepars;i++) delete [] curvature[i];
  delete [] curvature;
  delete(matrixcalc);
}

void CCF2SFit::CalcErrorMatrixFromCurvature(double **C){
  int i,j,k;
  CGSLMatrix_Real *matrixcalc;
  double **U;
  double *EigenVal;
  matrixcalc=new CGSLMatrix_Real(nfreepars);
  EigenVal=new double[nfreepars];
  U=new double *[nfreepars];
  for(i=0;i<nfreepars;i++) U[i]=new double[nfreepars];
  matrixcalc->EigenFind(C,U,EigenVal);
  
  for(i=0;i<npars;i++){
    for(j=0;j<=i;j++){
      if(i<nfreepars && j<nfreepars){
	ErrorMatrix[i][j]=0.0;
	for(k=0;k<nfreepars;k++){
	  ErrorMatrix[i][j]+=U[i][k]*U[j][k]/EigenVal[k];
	}
      }
      else ErrorMatrix[i][j]=0.0;
      if(i!=j) ErrorMatrix[j][i]=ErrorMatrix[i][j];
    }
    par[i]->error=sqrt(fabs(ErrorMatrix[i][i]));
  }

  for(i=0;i<nfreepars;i++) delete [] U[i];
  delete [] U;
  delete [] EigenVal;
  delete(matrixcalc);
}

#include "sfit_3dgaussian.cc"
#include "sfit_GX1d.cc"
#include "sfit_blast.cc"


#endif

