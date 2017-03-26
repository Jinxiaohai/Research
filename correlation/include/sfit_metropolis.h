#ifndef __INCLUDE_METROPOLIS_H__
#define __INCLUDE_METROPOLIS_H__

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string>

#include "coral.h"

class CParInfo{
 public:
  bool fixed;
  double x,xmin,xmax,error,oldx,xbar;
  char name[20];
  void Set(string sname,double xset,double error,
	   double xminset,double xmaxset);
};

class CCF2S_Metropolis{
public:
  CParInfo *par;
  int ndim;  // =1 or 3 for 1/3 dimensional CF
  int npars,nfreepars,Ncalls;
  double chisquared,oldchisquared;
  double *x_best;
  double chisquared_best;
  double **sigma;

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

  void Init();
  void CalcChiSquare();
  void UpdateStepSize();
  void Minimize(int maxcalls);
  void PrintPars();
  void PrintSigma();
  void PrintUsigma();
  void Reset();
  void FixPar(int ipar);
  void FreePar(int ipar);
  CCF2S_Metropolis();


 private:
  CRandom *randy;
  int Nsuccess;
  void UpdateSigma();
  CGSLMatrix_Real *matrix;
  double *xran;

 protected:
  double **Usigma;
  double *Ueigen;
};


class CCF2S_Metropolis_Blast : public CCF2S_Metropolis{
public:
  CCF2S_Metropolis_Blast(CSourceCalc *scset,C3DArray *cexpset,
		     C3DArray *cerrorset,C3DArray *ctheory3Dset,
		     CCHArray *ctheoryset,CCHArray *sourceset,
		     CKernel *kernelset);
};

class CCF2S_Metropolis_GX1D : public CCF2S_Metropolis{
 public:
  CCF2S_Metropolis_GX1D(CSourceCalc *scset,CCHArray *cexpset,
	       CCHArray *cerrorset,CCHArray *ctheoryset,
	       CCHArray *sourceset,CKernel *kernelset);
};

class CCF2S_Metropolis_3DGaussian : public CCF2S_Metropolis{
public:
  CCF2S_Metropolis_3DGaussian(CSourceCalc *scset,C3DArray *cexpset,
		     C3DArray *cerrorset,C3DArray *ctheory3Dset,
		     CCHArray *ctheoryset,CCHArray *sourceset,
		     CKernel *kernelset);
};


#endif
