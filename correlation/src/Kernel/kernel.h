#ifndef __INCLUDE_KERNEL_H
#define __INCLUDE_KERNEL_H

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <cstring>
using namespace std;

#include "sf.h"
#include "wavefunction.h"
#include "parametermap.h"
#include "misc.h"

class CKernel{
 public:
  double GetValue(int ell,double q, double r);
  double GetValue(int ell,int iq,int ir);
  void Read(char *datadirname);
  void Write(char *datadirname);
  void Print();
  int GetLMAX();
  double GetDELR();
  double GetDELQ();
  int GetNQMAX();
  int GetNRMAX();
  void Calc(CWaveFunction *wf);
  void Calc_ClassCoul(double ma,double mb,int zazb);
  void Calc_PureHBT();
  bool GetIDENTICAL();
  double GetPsiSquared(int iq,int ir,double ctheta);
  double GetPsiSquared(int iq,double r,double ctheta);
  double GetPsiSquared(double q,double r,double ctheta);
  CKernel(char *kparsfilename);
  ~CKernel();
 private:
  bool IDENTICAL;
  int ellmax;
  int nrmax,nqmax;
  double delr,delq;
  double ***kernel;
  double *P;
  void ParsInit(char *kparsfilename);
  double CalcPsiSquared(int iq,int ir,double ctheta);
  double CalcPsiSquared(int iq,double r,double ctheta);
  void CalcP(double ctheta);
};

class CKernelWF{
 public:
  double GetPsiSquared(int iq,int ir,int ictheta);
  double GetPsiSquared(int iq,int ir,double ctheta);
  double GetPsiSquared(int iq,double r,double ctheta);
  double GetPsiSquared(double q,double r,double ctheta);
  void Calc(CWaveFunction *wf);
  void Read(char *datadirname);
  void Write(char *datadirname);
  double GetDELR();
  double GetDELQ();
  double GetDELCTHETA();
  int GetNQMAX();
  int GetNRMAX();
  int GetNCTHETA();
  bool GetIDENTICAL();
  CKernelWF(char *kparsfilename);
  ~CKernelWF();

 private:
  bool IDENTICAL;
  int nctheta,nrmax,nqmax;
  double delr,delq,delctheta;
  double ***wfarray;
  void ParsInit(char *kparsfilename);
  
};

#endif
