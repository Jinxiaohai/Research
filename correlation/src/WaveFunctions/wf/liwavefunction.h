#ifndef __CWAVEFUNCTION_WF_H
#define __CWAVEFUNCTION_WF_H

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <cstring>

#include "misc.h"
#include "sf.h"
#include "parametermap.h"

using namespace std;

class CPlaneWave;
class CPartWave;

class CWaveFunction{
 public:
  int GetNQMAX();
  double GetQ(int iq);
  double GetDELQ();
  void PrintCdelta(double Rx,double Ry,double Rz);
  double GetPsiSquared(double *pa,double *xa,double *pb,double *xb);
  double GetPsiSquared(double q,double r,double ctheta);
  virtual double CalcPsiSquared(int iq,double r,double ctheta);
  CWaveFunction();
  ~CWaveFunction();
  double **Wepsilon,**delta,**ddeltadq,*eta,*channelweight;
  double GetIW(int ell,double epsilon,double q,int q1q2,
         double eta,double delta);
  CPlaneWave **planewave;
  CPartWave ***partwave;
 protected:
  void ParsInit(char *parsfilename);
  double PI;
  complex<double> ci;
  double HBARC;
  double MPI,MKAON,MPROTON,MLAMBDA,MNEUTRON;
  int *ell;
  //double **Wepsilon,**delta,**ddeltadq,*eta,*channelweight;
  int nchannels,q1q2;
  double m1,m2;
  bool generic;
  double mu,muscale,symmweight,epsilon;
  int q1q2scale;
  bool STRONG,COULOMB;
  //CPlaneWave **planewave;
  //CPartWave ***partwave;
  // These last two are redefined (overwritten) in inherited classes
  void InitArrays();
  void InitWaves();
  void getqrctheta(double *pa,double *xa,double *pb,double *xb,
		   double *q,double *r,double *ctheta);
  void EffectiveRange(int ichannel,double scattlength,double Reff);
  //double GetIW(int ell,double epsilon,double q,int q1q2,
  //     double eta,double delta);
  void phaseshift_CoulombCorrect(int ell,double q,double eta,
					      double &delta,double &ddeltadq);
  int ellmax,nqmax;
  double delq;
  double *qarray;
};

class CPlaneWave{
 public:
  CPlaneWave(double etaset,int Q1Q2,double qset);
  complex<double> planewave(double r,double ctheta);
  complex<double> hyper(complex<double> a,complex<double> b,
			complex<double> cz);
 private:
  double PI;
  complex<double>ci;
  double delx;
  int q1q2,nxmax;
  complex<double> couly;
  complex<double> cfact1;
  complex<double> cfact2;
  complex<double> chype[41];
  complex<double> chype1[6];
  complex<double> chype2[6];
  double q,eta; // q is the reduced mom., (p1-p2)/2
  complex<double> *hyperarray;
};

class CPartWave{
 public:
  CPartWave(double etaset,int q1q2,double qset,int ell,double epsilonset);
  CPartWave::CPartWave(double etaset,int q1q2set,double qset,int ellset,
		       double epsilonset,int nrmaxset,double delrset);
  ~CPartWave();
  complex<double> GetPhiIncoming(double r);
  complex<double> GetPhiOutgoing(double r);
  double sigma;
 private:
  complex<double> ci;
  double epsilon;
  double PI;
  int ell,q1q2;
  double q,eta;
  complex<double> *phi; /* An incoming Coulomb partial wave
			   which behaves as e^{-i(kr-eta*ln(2kr)+sigma)} */ 
  int nrmax;
  double delr;
  void phi_init();
};

// Parent class is CWaveFunction

class CWaveFunction_pkplus : public CWaveFunction {
 public:
  double CalcPsiSquared(int iq,double r,double ctheta);
  void read_phaseshifts();
  CWaveFunction_pkplus(char *parsfilename);
};
// Parent class is CWaveFunction

class CWaveFunction_ppiplus : public CWaveFunction {
 public:
  double CalcPsiSquared(int iq,double r,double ctheta);
  void read_phaseshifts();
  CWaveFunction_ppiplus(char *parsfilename);
};

class CWaveFunction_generic : public CWaveFunction {
 public:
  double CalcPsiSquared(int iq,double r,double ctheta);
  void reset(int q1q2,double m1,double m2,double symmweight);
  CWaveFunction_generic(char *parsfilename,int q1q2,double m1,double m2,double symmweight);
  
};

class CWaveFunction_Xipi : public CWaveFunction {
 public:
  double CalcPsiSquared(int iq,double r,double ctheta);
  void get_phaseshifts_Xistar(double q,double &delta,double &ddeltadq);
  CWaveFunction_Xipi(char *parsfilename);
};

class CWaveFunction_pn : public CWaveFunction {
 public:
  double CalcPsiSquared(int iq,double r,double ctheta);
  //void get_phaseshifts_pn(double q,double &delta,double &ddeltadq);
  void read_phaseshifts();
  CWaveFunction_pn(char *parsfilename);
};

class CWaveFunction_plambda : public CWaveFunction {
 public:
  double CalcPsiSquared(int iq,double r,double ctheta);
  //void get_phaseshifts_plambda(double q,double &delta,double &ddeltadq);
  void get_phaseshifts();
  CWaveFunction_plambda(char *parsfilename);
};

class CWaveFunction_kpluspiminus : public CWaveFunction {
 public:
  double CalcPsiSquared(int iq,double r,double ctheta);
  void get_phaseshifts_kpluspiminus();
  CWaveFunction_kpluspiminus(char *parsfilename);
};

class CWaveFunction_pp : public CWaveFunction {
 public:
  double CalcPsiSquared(int iq,double r,double ctheta);
  //void get_phaseshifts_pp(double q,double &delta,double &ddeltadq);
  void read_phaseshifts();
  CWaveFunction_pp(char *parsfilename);
};

class CWaveFunction_nn : public CWaveFunction {
 public:
  double CalcPsiSquared(int iq,double r,double ctheta);
  //void get_phaseshifts_nn(double q,double &delta,double &ddeltadq);
  void read_phaseshifts();
  CWaveFunction_nn(char *parsfilename);
};

class CWaveFunction_pipluspiminus : public CWaveFunction {
 public:
  double CalcPsiSquared(int iq,double r,double ctheta);
  CWaveFunction_pipluspiminus(char *parsfilename);
};

class CWaveFunction_pipluspiplus : public CWaveFunction {
 public:
  double CalcPsiSquared(int iq,double r,double ctheta);
  CWaveFunction_pipluspiplus(char *parsfilename);
};

class CWaveFunction_lambdalambda : public CWaveFunction {
 public:
  double CalcPsiSquared(int iq,double r,double ctheta);
  CWaveFunction_lambdalambda(char *parsfilename);
  void GetPhaseshifts();
};

class CWaveFunction_lambdalambda_antiparspin : public CWaveFunction {
 public:
  double CalcPsiSquared(int iq,double r,double ctheta);
  CWaveFunction_lambdalambda_antiparspin(char *parsfilename);
  void GetPhaseshifts();
};

class CWaveFunction_lambdalambda_parspin : public CWaveFunction {
 public:
  double CalcPsiSquared(int iq,double r,double ctheta);
  CWaveFunction_lambdalambda_parspin(char *parsfilename);
  void GetPhaseshifts();
};

class CWaveFunction_kpi : public CWaveFunction {
public:
  double CalcPsiSquared(int iq,double r,double ctheta);
  CWaveFunction_kpi(char *parsfilename);
};

class	CWaveFunction_pps : public CWaveFunction
{
 public:
  double CalcPsiSquared(int iq,double r,double ctheta);
  CWaveFunction_pps(char *parsfilename);

	bool	bReid;
	double	fRMax;
	long	nRMax;

	~CWaveFunction_pps();

	char	pname[4];
	char	ptype[3];
	short	nl;
	short	nj;

	complex <double>	*	*	pcDelPsi10;	// DelPsi(|1,0>) = DelPsi(|-1,0>)
	complex <double>	*	*	pcDelPsi01;	// DelPsi(|0,1>) = DelPsi(|0,-1>)
	complex <double>	*	*	pcDelPsi11;	// DelPsi(|1,-1>) = DelPsi(|-1,1>)
	complex <double>	*	*	pcDelPsi00;	// DelPsi(|0,0>)
	complex <double>	*	*	pcDelPsiSinglet; // Singlet State

	double	*	pfPhaseShift[4];

	double	*	pfPhaseShiftNS;		//No-Strong Interaction PhaseShift
	
	double	*	pfPhaseShift3P2NoF;	// no coupling to F state

	bool	bCoupleToF;

protected:
	bool Initialize(void);

	bool GetPsiNoneCoupledState(long iq, short n_L, const char * szPartialWaveName, short nTagPartialWave, bool bReidInput);

	bool GetPsiCoupledStates(long iq, const char * szPartialWaveName);

	complex <double>	*	pcU1;
	complex <double>	*	pcU2;
	complex <double>	*	pcW1;
	complex <double>	*	pcW2;
	complex <double>	*	pcPsi;

	complex <double>	*	pc3P0;	// None Coupled State
	complex <double>	*	pc3P1;	// None Coupled State

	complex <double>	*	pc3P2;	// Coupled State with L=J-1
	complex <double>	*	pc3F2;	// Coupled State with L=J+1

	complex <double>	*	pc3P2S;	// Spin-flipped term
	complex <double>	*	pc3F2S;	// Spin-flipped term

	complex <double>	*	pc3P2NoF;	// No coupling to F state
    

	complex <double>	*	pc1S0;	// Singlet State
	complex <double>	*	pcNS;	//No-Strong Interaction Wavefunctions

	static	const	short	TAG_3P0 = 1;
	static	const	short	TAG_3P1 = 2;
	static	const	short	TAG_1S0 = 3;
	static	const	short	TAG_NS = 4;
	static	const	short	TAG_3P2NOF = 5;

};


#endif
