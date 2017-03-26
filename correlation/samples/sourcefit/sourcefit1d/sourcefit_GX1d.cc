#include "coral.h"

int main (int argc, char *argv[]){
  double lambdaG=0.4,R=4,lambdaX=0.2,X=10.0,a=5.0;
  int lx,ly,lz;
  double *xcontour,*ycontour;
  int nfound_contour=0,npts_contour=20;
  xcontour=new double[npts_contour];
  ycontour=new double[npts_contour];

  CCF2SFit_GX1D *fitter;
  CCHArray *Ctheory;
  CCHArray *S;
  CCHArray *Cdata;
  CCHArray *Cerror;
  CSourceCalc_GX1D *scalc;
  CKernel *kernel;
  char filename[200];
  char dirname[200];
  double chisquare;
  char pardirname[100];
  sprintf(pardirname,"parameters/pp\0");

  sprintf(filename,"%s/aparsCH_source.dat\0",pardirname);
  S=new CCHArray(filename);
  sprintf(filename,"%s/aparsCH_CF.dat\0",pardirname);
  Ctheory=new CCHArray(filename);
  Cdata=new CCHArray(filename);
  sprintf(dirname,"fakedata/GX1d/pp\0");
  Cdata->ReadAX(dirname);
  Cerror=new CCHArray(filename);
  sprintf(dirname,"fakedata/GX1d/pp_errors\0");
  Cerror->ReadAX(dirname);

  sprintf(filename,"%s/kparameters.dat\0",pardirname);
  kernel=new CKernel(filename);
  sprintf(dirname,"../kdata/pp\0");
  kernel->Read(dirname);
  
  scalc=new CSourceCalc_GX1D();
  scalc->SetSPars(0.5,6.0,0.0,10.0,5.0);
  scalc->CalcS(S);
  scalc->NormCheck(S);

  fitter=new CCF2SFit_GX1D(scalc,Cdata,Cerror,Ctheory,S,kernel);
  lx=ly=lz=0;
  fitter->SetL(lx,ly,lz);
  fitter->FixPar("lambdaX");
  fitter->FixPar("a");
  fitter->FixPar("X");

  // Find minimum
  printf("Beginning Parameters : \n");
  fitter->PrintPars();
  fitter->Newton(6);
  fitter->PrintPars();

  return 0;
}
