#ifndef __INCLUDE_SOURCECALC_H
#define __INCLUDE_SOURCECALC_H
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <fstream>
#include "arrays.h"
#include "parametermap.h"
#include "random.h"

using namespace std;

class CSourceCalc{
 public:
  parametermap spars;
  virtual void CalcS(CCHArray *A);
  virtual void CalcS(int lx,int ly,int lz,CCHArray *A);
  void ReadSPars(char *sparsfilename);
  void NormCheck(CCHArray *A);
  void NormCheck(C3DArray *threed);
  void CalcEffGaussPars(CCHArray *A);
  void CalcEffGaussPars(CCHArray *A,double &Rx,double &Ry,double &Rz,
			double &Xoff,double &Yoff,double &Zoff);
  CSourceCalc::CSourceCalc();
};

class CSourceCalc_GX1D : public CSourceCalc{
/*
  S = { N_G * lambda_G * exp(-r^2/4R^2)
        + N_X * lambda_X * exp(-[r^2/X^2+4R^4/X^4]) } * [r^2/(r^2+a^2)]^(L/2)
	where N_x and N_g are constants that make individual contributions
        integrate to unity for L=0
*/
 public:
  CSourceCalc_GX1D();
  void InitSPars();
  void SetSPars(double lambdaG,double R,double lambdaX,double X,double a);
  void CalcS(int lx,int ly,int lz,CCHArray *A);
  void CalcS(CCHArray *A);
};

class CSourceCalc_Gaussian : public CSourceCalc{
 public:
  CSourceCalc_Gaussian();
  void SetSPars(double lambdaset,double Rxset,double Ryset,double Rzset);
  void SetSPars(double lambdaset,double Rxset,double Ryset,double Rzset,
		double Xoffset,double Yoffset,double Zoffset);
  void SetSPars(double lambdaset,double Rxset,double Ryset,double Rzset,
		double Xoffset,double Yoffset,double Zoffset,
		double Euler_phiset,double Euler_thetaset,
		double Euler_psiset);
  void CalcS(CCHArray *A);
  void CalcS(C3DArray *threed);
  
 private:
  void CalcAlpha(double **alpha,CCHArray *A);
  void InitSPars();
  // S ~ exp{-alpha_ij (x_i-off_i)(x_j-off_j)} 
  // Calcs ignore zeroth components
};

class CSourceCalc_EllipticBlast : public CSourceCalc{
 public:
  CSourceCalc_EllipticBlast();
  void CalcS(CCHArray *A);
  void CSourceCalc_EllipticBlast::SetSPars(double Rxset,double Ryset,
					   double Tauset,
					   double BetaXset,double BetaYset,
					   double Tset,double Ptset,
					   double Phiset,double EtaGset,
					   double Maset,double Mbset);
  void CSourceCalc_EllipticBlast::SetSPars(double Rset,double Tauset,
			       double Betaset,double Tset,double Ptset);
  CRandom *randy;
 private:
  void CSourceCalc_EllipticBlast::Get_r(double *p,int nsample,double **r);
  void InitSPars();
};

class CSourceCalc_Blast : public CSourceCalc{
 public:
  CSourceCalc_Blast();
  void CalcS(CCHArray *A);
  void CSourceCalc_Blast::SetSPars(double lambdaset,
				   double Rset,double Tauset,double DelTauset,
				   double Betaset,double Tset,double Ptset,
				   double EtaGset,double Maset,double Mbset);
  void CSourceCalc_Blast::SetSPars(double lambdaset,
				   double Rset,double Tauset,double DelTauset);
  CRandom *randy;
 private:
  void CSourceCalc_Blast::Get_r(double *p,double **r);
  void InitSPars();
  double GetTau(double tau0,double deltau);
};

class CSourceCalc_OSCAR : public CSourceCalc{
 public:
  CSourceCalc_OSCAR();
  void CalcS(CCHArray *A);
  void SetSPars(double Pt_set,double delpt_set,
		double phimin_deg_set,double phimax_deg_set,
		double ymin_set,double ymax_set);
  void CSourceCalc_OSCAR::SetIDs(int IDa_set,int IDb_set);
  
 private:
  void ReadR(double **ra,double **rb,int &na,int &nb);
  void InitSPars();
};

#endif
