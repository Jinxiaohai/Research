#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iostream>

using namespace std;

#include "coral.h"

int main(){
  CCHArray *Cfake;
  C3DArray *Cfake3D;
  C3DArray *Ctheory3D;
  C3DArray *Efake3D;
  CCHArray *Sgauss;
  CSourceCalc_Gaussian *scalc;
  CWaveFunction *wf;
  CKernel *kernel;
  double lambda=0.5,Rx=7.0,Ry=5.0,Rz=3.0,Xoff=0.0;
  double error=0.01,chisquare;
  char pardirname[100];
  sprintf(pardirname,"../parameters/pp\0");

  // Create Arrays
  char apars_filename[120];
  sprintf(apars_filename,"%s/aparsCH_source.dat\0",pardirname);
  Sgauss=new CCHArray(apars_filename);
  sprintf(apars_filename,"%s/aparsCH_CF.dat\0",pardirname);
  Cfake=new CCHArray(apars_filename);
  sprintf(apars_filename,"%s/apars3D_CF.dat\0",pardirname);
  Cfake3D=new C3DArray(apars_filename);
  Ctheory3D=new C3DArray(apars_filename);
  Efake3D=new C3DArray(apars_filename);

  char fakedatadirname[120];
  sprintf(fakedatadirname,"3dgaussian/pp\0");
  char fakeerrordirname[120];
  sprintf(fakeerrordirname,"3dgaussian/pp_errors\0");

  char wfparsfilename[120];
  sprintf(wfparsfilename,"%s/kparameters.dat\0",pardirname);
  char kdatadirname[200];
  sprintf(kdatadirname,"../../kdata/pp\0");

  kernel=new CKernel(wfparsfilename);
  kernel->Read(kdatadirname);
  //wf=new CWaveFunction_pp(wfparsfilename);
  //kernel->Calc(wf);
  //kernel->Write(kdatadirname);

  scalc=new CSourceCalc_Gaussian();
  scalc->SetSPars(lambda,Rx,Ry,Rz,Xoff,0.0,0.0);
  scalc->CalcS(Sgauss);
  scalc->NormCheck(Sgauss);
  S2CF::s2c(Sgauss,kernel,Cfake);
  ArrayCalc::Calc3DArrayFromAExpArray(Cfake,Ctheory3D);

  //Cfake3D->ReadArray(fakedatadirname);
  //Efake3D->ReadArray(fakeerrordirname);

  // Add random errors to Cfake3D
  Efake3D->Randomize_Gaussian();
  Efake3D->ScaleArray(error);
  ArrayCalc::AddArrays(Ctheory3D,Efake3D,Cfake3D);
  // Make error matrix and calculate chi square
  Efake3D->MakeConstant(error);
  chisquare=CFCalc::GetChiSquared(Cfake3D,Efake3D,Ctheory3D);
  printf("chisquare=%g\n",chisquare);
  
  // Write Arrays
  Cfake3D->WriteArray(fakedatadirname);
  Efake3D->WriteArray(fakeerrordirname);

  return 0;
  
}

