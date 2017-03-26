#ifndef __INCLUDE_ARRAYS_H
#define __INCLUDE_ARRAYS_H

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include "random.h"
#include "parametermap.h"
#include "sf.h"

using namespace std;

class CCHArray{
 public:
  CCHArray(char *arrayparsfilename);
  CCHArray(int LMAXset,int NRADIALset,double RADSTEPset);
  CCHArray(int LMAXset,int NRADIALset,double RADSTEPset,
		 bool XSYMset,bool YSYMset,bool ZSYMset);
  ~CCHArray();

  int GetLMAX();
  int GetNRADIAL();
  double GetRADSTEP();
  void SetLMAX(int LMAXset);
  void SetRADSTEP(double RADSTEPset);
  double GetElement(int lx,int ly,int lz,int ir);
  double GetElement(int lx,int ly,int lz,double r);
  void SetElement(int lx,int ly,int lz,int ir,double Element);
  void SetElement(int lx,int ly,int lz,double r,double Element);
  void IncrementElement(int lx,int ly,int lz,int ir,double increment);
  void IncrementElement(int lx,int ly,int lz,double r,double increment);
  void ScaleArray(double scalefactor);
  void ScaleArray(double scalefactor,int ir);
  void ZeroArray();
  void ZeroArray(int lx,int ly,int lz);
  void ZeroArray(int ir);
  void ZeroArray_Partial(int LMAX_Partial);
  void ZeroArray_Partial(int LMAX_Partial,int ir);
  void PrintArrayFixedIR(int ir);
  void PrintArrayFixedIR(int LMAXPrint,int ir);
  //FillRemainder functions use identity
  //   A_{lx+2,ly,lz}+A_{lx,ly+2,lz}+A_{lx,ly,lz+2}=0
  // to find entire array if (2L+1) values for lx=0,1 are known
  double GetBiggest(int ir);
  bool GetXSYM();
  bool GetYSYM();
  bool GetZSYM();
  void ReadAX(char *dirname);
  void WriteAX(char *dirname);
  void ReadAllA(char *dirname);
  void WriteAllA(char *dirname);
  void WriteShort(char *filename,int WLMAX);
  void PrintPars();
  void IncrementAExpArray(double x,double y,double z,double weight);
  void IncrementAExpArrayFromE(double ex,double ey,double ez,
			       double weight,int ir);
  void AltIncrementAExpArrayFromE(double ex,double ey,double ez,
				  double weight,int ir);
  void AltAltIncrementAExpArrayFromE(double ex,double ey,double ez,
				  double weight,int ir);
  void IncrementMArrayFromE(double ex,double ey,double ez,double weight,int ir);
  void IncrementAExpArrayFromThetaPhi(double theta,double phi,
				      double weight,int ir);
  void IncrementMArrayFromThetaPhi(double theta,double phi,
				   double weight,int ir);
  double GetMElementFromAExpArray(int lx,int ly,int lz,int ir);
  double GetAExpElementFromMArray(int lx,int ly,int lz,int ir);
  void FillRemainderX(int ir);
  void FillRemainderY(int ir);
  void FillRemainderZ(int ir);
  void FillRemainderX();
  void FillRemainderY();
  void FillRemainderZ();
  double AExpand(double ex,double ey,double ez,int ir);
  double AExpand(double x,double y,double z);
  double AExpand(double theta,double phi,int ir);
  void Detrace(int ir);
  void Detrace();
  void RandomInit(int iseed);
  void Randomize(double mag,int ir);
  void RandomizeA(double mag,int ir);  
  void Randomize(double mag);
  void RandomizeA(double mag);  
  void RandomizeA_Gaussian(double mag,int ir);
  void RandomizeA_Gaussian(double mag);

 private:
  static CCHCalc *chcalc;
  static CRandom *randy;
  bool XSYM,YSYM,ZSYM;
  double RADSTEP;
  int NRADIAL;
  int LMAX;
  double ****A;
  int dlx,dly,dlz;
  void CreateArray();

};


class C3DArray{
 public:
  C3DArray(char *arrayparsfilename);
  C3DArray(int NXYZMAX,double DELXYZ,bool XSYM,bool YSYM,bool ZSYM);
  C3DArray(int NXMAX,double DELX,int NYMAX,double DELY,int NZMAX,double DELZ,
	   bool XSYM,bool YSYM,bool ZSYM);
  ~C3DArray();
  int GetNXMAX();
  int GetNYMAX();
  int GetNZMAX();
  double GetDELX();
  double GetDELY();
  double GetDELZ();
  double GetElement(double x,double y,double z);
  double GetElement_Interpolate(double x,double y,double z);
  double GetElement(int isx,int ix,int isy,int iy,int isz,int iz);
  void SetElement(int isx,int ix,int isy,int iy,int isz,int iz,double value);
  void IncrementElement(int isx,int ix,int isy,int iy,int isz,int iz,
			  double value);
  void SetElement(double x,double y,double z,double value);
  void IncrementElement(double x,double y,double z,double increment);
  void ScaleArray(double scalefactor);
  void ZeroArray();
  void PrintArray();
  double GetBiggest();
  bool GetXSYM();
  bool GetYSYM();
  bool GetZSYM();
  void ReadArray(char *dirname);
  void WriteArray(char *dirname);
  void PrintPars();
  void Randomize();
  void Randomize_Gaussian();
  void MakeConstant(double c);
 private:
  bool XSYM,YSYM,ZSYM;
  double DELX,DELY,DELZ;
  int NXMAX,NYMAX,NZMAX;
  double ******F;
  void ReadPars(char *arrayparsfilename);
  void CreateArray();
  void DeleteArray();
  static CRandom *randy;
};

class CYlmArray{
 public:
  CYlmArray(int LMAXset,int NRADIALset);
  ~CYlmArray();
  int GetLMAX();
  complex<double> GetElement(int L,int M,int ir);
  void SetElement(int L,int M,int ir,complex<double>);
  void IncrementElement(int L,int M,int ir,complex<double> increment);
  void ScaleArray(double scalefactor);
  void ScaleArray(double scalefactor,int ir);
  void ZeroArray();
  void ZeroArray(int ir);
  void PrintArrayFixedIR(int ir);
  void PrintArrayFixedIR(int LMAXPrint,int ir);
 private:
  int NRADIAL;
  int LMAX;
  complex<double> ***ylm;
};

namespace ArrayCalc{

  void CalcMArrayFromAExpArray(CCHArray *A,CCHArray *M);
  void CalcMArrayFromAExpArray(CCHArray *A,int ira,CCHArray *M,int irm);
  void CalcAExpArrayFromMArray(CCHArray *M,CCHArray *A);
  void CalcAExpArrayFromMArray(CCHArray *M,int irm,CCHArray *A,int ira);

  void AddArrays(CCHArray *A,CCHArray *B,CCHArray *C);
  void AddArrays(CCHArray *A,int ira,CCHArray *B,int irb,CCHArray *C,int irc);
  // If A(Omega)=B(Omega)*C(Omega), this finds C in terms of A and B
  void DivideArrays(CCHArray *A,CCHArray *B,CCHArray *C);
  void DivideArrays(CCHArray *A,int ira,CCHArray *B,int irb,
		    CCHArray *C,int irc);
  // If C(Omega)=A(Omega)*B(Omega), this finds A in terms of A and B
  void MultiplyArrays(CCHArray *A,CCHArray *B,CCHArray *C);
  void MultiplyArrays(CCHArray *A,int ira,CCHArray *B,
		      int irb,CCHArray *C,int irc);
  // If you know array is zero up to given Ls, or don't care to detrace,
  // this can save time
  void MultiplyArrays_Partial(int LMAXA,CCHArray *A,int ira,
			      int LMAXB,CCHArray *B,int irb,
			      int LMAXC,CCHArray *C,int irc);

  void CopyArray(CCHArray *A,CCHArray *B);
  void CopyArray(CCHArray *A,int ira,CCHArray *B,int irb);

  void CalcYlmExpArrayFromAExpArray(CCHArray *A,int ir,
				    CYlmArray *YlmArray,int irlm);
  void CalcAExpArrayFromYlmExpArray(CYlmArray *YlmArray,int irlm,
				    CCHArray *A,int ira);

  // If one has an angular function
  // F=\sum_{\vec\ell} M_{\vec\ell} e_x^{\ell_x}e_y^{\ell_y}e_z^{\ell_z}
  // and wants to find expansion coefficients A which give the same answer
  // but satisfy traceless condition, this finds the array A
  void Detrace(CCHArray *M,CCHArray *A);
  void Detrace(CCHArray *M,int irm,CCHArray *A,int ira);

  // XExpArray express functions in form exp( X_ell * nhat^ell ) 
  void CalcAExpArrayFromXExpArray(CCHArray *X,CCHArray *A);
  void CalcAExpArrayFromXExpArray(CCHArray *X,int irx,CCHArray *A,int ira);
  void CalcXExpArrayFromAExpArray(CCHArray *A,CCHArray *X);
  void CalcXExpArrayFromAExpArray(CCHArray *A,int ira,CCHArray *X,int irx);


  void CalcAExpArrayFrom3DArray(C3DArray *threedarray,CCHArray *A);
  void Calc3DArrayFromAExpArray(CCHArray *A,C3DArray *threed);

  void MultiplyArrays(C3DArray *A,C3DArray *B,C3DArray *C);
  void DivideArrays(C3DArray *A,C3DArray *B,C3DArray *C);
  void AddArrays(C3DArray *A,C3DArray *B,C3DArray *C);
  void CopyArray(C3DArray *A,C3DArray *B);

  void InvertArray(C3DArray *A,C3DArray *B);
  void InvertArray(CCHArray *A,int ira,CCHArray *B,int irb);
  void InvertArray(CCHArray *A,CCHArray *B);

  /* These check that two arrays have the same parameters
     For arrays of different types (e.g., 3D and Cart.Harmonic) it checks that
     they have the same reflection symmetries */
  bool CompareArrayParameters(C3DArray *threed,CCHArray *A);
  bool CompareArrayParameters(CCHArray *A,C3DArray *threed);
  bool CompareArrayParameters(CCHArray *A,CCHArray *B);
  bool CompareArrayParameters(C3DArray *threeda,C3DArray *threedb);
};

#endif

