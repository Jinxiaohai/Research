#include "coral.h"

int main (int argc, char *argv[]){
  double *xcontour,*ycontour;
  int nfound_contour=0,npts_contour=20;
  xcontour=new double[npts_contour];
  ycontour=new double[npts_contour];

  CCF2SFit *fitter;
  CCHArray *Ctheory;
  CCHArray *Sgauss;
  C3DArray *Ctheory3D;
  C3DArray *Cdata3D;
  C3DArray *Cerror3D;
  CSourceCalc_Gaussian *scalc;
  CKernel *kernel;
  char filename[200];
  char dirname[200];
  double chisquare;
  char pardirname[100];
  sprintf(pardirname,"parameters/pp\0");

  sprintf(filename,"%s/aparsCH_source.dat\0",pardirname);
  Sgauss=new CCHArray(filename);
  sprintf(filename,"%s/aparsCH_CF.dat\0",pardirname);
  Ctheory=new CCHArray(filename);
  sprintf(filename,"%s/apars3D_CF.dat\0",pardirname);
  Ctheory3D=new C3DArray(filename);
  
  Cdata3D=new C3DArray(filename);
  sprintf(dirname,"fakedata/3dgaussian/pp\0");
  Cdata3D->ReadArray(dirname);
  Cerror3D=new C3DArray(filename);
  sprintf(dirname,"fakedata/3dgaussian/pp_errors\0");
  Cerror3D->ReadArray(dirname);

  sprintf(filename,"%s/kparameters.dat\0",pardirname);
  kernel=new CKernel(filename);
  sprintf(dirname,"../kdata/pp\0");
  kernel->Read(dirname);

  scalc=new CSourceCalc_Gaussian();
  scalc->SetSPars(0.5,7,5,3,0,0,0);
  scalc->CalcS(Sgauss);
  scalc->NormCheck(Sgauss);

  fitter=new CCF2SFit_3DGaussian(scalc,Cdata3D,Cerror3D,Ctheory3D,
				 Ctheory,Sgauss,kernel);
  fitter->FixPar("Xoff");

  // Find minimum
  fitter->PrintPars();

  //fitter->Newton(6);
  fitter->Metropolis(50);
  fitter->UpdateMetroStepSize();
  fitter->Metropolis(50);


  fitter->PrintPars();

  for(int iq=0;iq<Ctheory->GetNRADIAL();iq++)
    printf("C(%d)=%g\n",iq,Ctheory->GetElement(0,0,0,iq));

  return 0;
}
