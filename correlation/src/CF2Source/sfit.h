#ifndef __INCLUDE_SFIT_H__
#define __INCLUDE_SFIT_H__

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string>

#include "coral.h"

class CParInfo{
 public:
  bool fixed;
  double xmin,xmax,error,xbar;
  double bestx,currentx;
  char *name;
  void Set(string sname,double xset,double error,
	   double xminset,double xmaxset);
  CParInfo();
  ~CParInfo();
};


class CCF2SFit{
 public:
  void SetPar(string parstring,double value);
  void PrintPars();
  void PrintErrorMatrix();
  void PrintStepMatrix();
  void FixPar(string parname);
  void FreePar(string parname);
  void UseBestPars();
  void SetL(int lxset,int lyset,int lzset);

  void Metropolis(int maxcalls);
  void SteepestDescent(int maxtries);
  void Newton(int maxtries);
  void UpdateMetroStepSize();

  CCF2SFit();
  ~CCF2SFit();

 protected:
  CParInfo **par;
  int ndim;  // =1 or 3 for 1/3 dimensional CF
  int npars,nfreepars;
  double **ErrorMatrix;
  double currentchisquared,bestchisquared;

  CSourceCalc *sourcecalc;
  CKernel *kernel;

  CCHArray *ctheory;
  // These are used for 3D Calculations
  C3DArray *cexp3D;
  C3DArray *cerror3D;
  C3DArray *ctheory3D;
  //These are used for 1D Calculations
  int lx,ly,lz;
  CCHArray *cexp;
  CCHArray *cerror;
  CCHArray *source;

  CRandom *randy;
  void SwitchPars(int ipara,int iparb);

  int ncalls;
  double **StepMatrix;
  void Init();
  double GetChiSquared(double *x);
  void CalcErrorMatrixFromCurvature(double **C);
};


class CCF2SFit_Blast : public CCF2SFit{
 public:
  CCF2SFit_Blast(CSourceCalc *scset,C3DArray *cexpset,
		 C3DArray *cerrorset,C3DArray *ctheory3Dset,
		 CCHArray *ctheoryset,CCHArray *sourceset,
		 CKernel *kernelset);
};

class CCF2SFit_GX1D : public CCF2SFit{
 public:
  CCF2SFit_GX1D(CSourceCalc *scset,CCHArray *cexpset,
		CCHArray *cerrorset,CCHArray *ctheoryset,
		CCHArray *sourceset,CKernel *kernelset);
};

class CCF2SFit_3DGaussian : public CCF2SFit{
 public:
  CCF2SFit_3DGaussian(CSourceCalc *scset,C3DArray *cexpset,
		      C3DArray *cerrorset,C3DArray *ctheory3Dset,
		      CCHArray *ctheoryset,CCHArray *sourceset,
		      CKernel *kernelset);
};


#endif

