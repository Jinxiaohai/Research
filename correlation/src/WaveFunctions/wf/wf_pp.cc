//CWaveFunction_pp.cpp
// Date: 20060621
// Here calculate l=1, P wave. Not tested for higher l yet.

#include "reid93.cc"
#include "wavefunction.h"
using	namespace	std;

CWaveFunction_pp::CWaveFunction_pp(char * parsfilename) : CWaveFunction()
{
  ParsInit(parsfilename);
  bReid = STRONG;
  nRMax = 1000;
  fRMax = 10.0;

  m1=MPROTON; 
  m2=m1; 
  q1q2=1; 
	
  nchannels=4;

  ellmax=1;
  InitArrays();
  printf("Arrays Initialized\n");
  ell[0]=0;
  ell[1]=ell[2]=ell[3]=1;

  InitWaves();
  printf("Partial Waves Initialized\n");

  //w.f. for PP scattering
  strcpy(ptype, "PP");

  nl = ellmax;

  pc3P0 = new complex<double> [nRMax+1];
  pc3P1 = new complex<double> [nRMax+1];
  pc1S0 = new complex<double> [nRMax+1];

  pc3P2NoF = new complex<double> [nRMax + 1];
  pfPhaseShift3P2NoF = new double[nqmax];

  pcNS = new complex<double> [nRMax+1];
    
  pcDelPsi10 = new complex<double> *[nqmax];
  pcDelPsi01 = new complex<double> *[nqmax];
  pcDelPsi11 = new complex<double> *[nqmax];
  pcDelPsi00 = new complex<double> *[nqmax];
  pcDelPsiSinglet = new complex<double> *[nqmax];

  for(int j=0;j<4;j++)
    {
      pfPhaseShift[j] = new double[nqmax];
    }

  pfPhaseShiftNS = new double[nqmax];


  for(long iq=0;iq<nqmax;iq++)
    {
      pcDelPsi10[iq] = new complex<double>[nRMax + 1];
      pcDelPsi01[iq] = new complex<double>[nRMax + 1];
      pcDelPsi11[iq] = new complex<double>[nRMax + 1];
      pcDelPsi00[iq] = new complex<double>[nRMax + 1];
      pcDelPsiSinglet[iq] = new complex<double>[nRMax + 1];
    }
  pcPsi = new complex<double>[nRMax + 1];

  Initialize();

}

bool CWaveFunction_pp::Initialize(void)
{
  for(long iq=0;iq<nqmax;iq++)
    {
      GetPsiNoneCoupledState(iq, 0, "1S0", TAG_1S0, bReid);
      GetPsiNoneCoupledState(iq, 1, "3P0", TAG_3P0, bReid);
      GetPsiNoneCoupledState(iq, 1, "3P1", TAG_3P1, bReid);
	  GetPsiNoneCoupledState(iq, 1, "3C2", TAG_3P2NOF, bReid);
      GetPsiNoneCoupledState(iq, 1, "NS", TAG_NS, false);

	  for(long ir=nRMax;ir>=0;ir--)
	    {
	      pcDelPsi10[iq][ir] = 0.5 * (pc3P2NoF[ir] - pc3P1[ir]);
	      pcDelPsi01[iq][ir] = 0.5 * (pc3P2NoF[ir] + pc3P1[ir] - 2.0 * pcNS[ir]);
	      pcDelPsi11[iq][ir] = 1.0 / 3.0 * (pc3P2NoF[ir] - pc3P0[ir]);
	      pcDelPsi00[iq][ir] = 1.0 / 3.0 * (2.0 * pc3P2NoF[ir] + pc3P0[ir] - 3.0 * pcNS[ir]);
	      pcDelPsiSinglet[iq][ir] = pc1S0[ir] - pcNS[ir];
	    }
    }
  return true;
}


CWaveFunction_pp::~CWaveFunction_pp()
{
	for(long iq=0;iq<nqmax;iq++)
	{
		delete[]	pcDelPsi10[iq];
		delete[]	pcDelPsi01[iq];
		delete[]	pcDelPsi11[iq];
		delete[]	pcDelPsi00[iq];
		delete[]	pcDelPsiSinglet[iq];

	}
	for(int j=0;j<4;j++)
	{
		delete[]	pfPhaseShift[j];
	}
	delete[]	pfPhaseShiftNS;
	delete[]	pcDelPsi10;
	delete[]	pcDelPsi01;
	delete[]	pcDelPsi11;
	delete[]	pcDelPsi00;
	delete[]	pcDelPsiSinglet;
	delete[]	pc3P0;
	delete[]	pc3P1;
	delete[]	pc1S0;
	delete[]	pc3P2NoF;
	delete[]	pfPhaseShift3P2NoF;
	delete[]	pcNS;

	delete[]	pcPsi;
}


double CWaveFunction_pp::CalcPsiSquared(int iq,double r,double ctheta)
{
  double psisquared;
  double q = (iq + 0.5) * delq;
  long ir = (long)(r / fRMax * nRMax);
  complex<double> psi,psisymm,psianti,psia,psib;
	
  psia=planewave[iq]->planewave(r,ctheta);
  psib=planewave[iq]->planewave(r,-ctheta);
  psisymm=(1.0/sqrt(2.0))*(psia+psib);
  psianti=(1.0/sqrt(2.0))*(psia-psib);

  if(STRONG==1)
    {
      if(r<epsilon)
	{
		//S and P waves
	  psisquared = 0.5 * real((psianti + pcDelPsi01[iq][ir]) * conj(psianti + pcDelPsi01[iq][ir])) + 0.5 * real(pcDelPsi10[iq][ir] * conj(pcDelPsi10[iq][ir]))
	    + 0.25 * real((psianti + pcDelPsi00[iq][ir]) * conj(psianti + pcDelPsi00[iq][ir])) + 0.25 * 2 * real(pcDelPsi11[iq][ir] * conj(pcDelPsi11[iq][ir]))
	    + 0.25 * real((psisymm + pcDelPsiSinglet[iq][ir]) * conj(psisymm + pcDelPsiSinglet[iq][ir]));

		//S wave only
	  //psisquared = 0.5 * real((psianti) * conj(psianti))
	  //  + 0.25 * real((psianti) * conj(psianti))
	  //  + 0.25 * real((psisymm + pcDelPsiSinglet[iq][ir]) * conj(psisymm + pcDelPsiSinglet[iq][ir]));
	}
      else
	{
	  double	theta=acos(ctheta);
	  double	x=q*r/HBARC;
	  double delta_3p0,delta_1s0,delta_3p1,delta_3p2;
	  complex<double> hstar0,hstar1;
	  complex<double> Xlm00,Xlm10,Xlm11;

	  delta_1s0=pfPhaseShift[0][iq] - 0.5 * PI;
	  delta_3p0=pfPhaseShift[1][iq] - 0.5 * PI;
	  delta_3p1=pfPhaseShift[2][iq] - 0.5 * PI;
      delta_3p2=pfPhaseShift3P2NoF[iq] - 0.5 * PI;

		// Only S wave
	  //delta_3p0=0.0;
	  //delta_3p1=0.0;
	  //delta_3p2=0.0;

	  hstar0=partwave[0][iq]->GetPhiIncoming(r)/x;
	  hstar1=partwave[1][iq]->GetPhiIncoming(r)/x;

	  // XlmLM 9s Y_{LM}*(1/2)*i^L*sqrt(4*PI*(2*L+1))*hstar_L
	  Xlm00=0.5*sqrt(2.0)*sqrt(4.0*PI)*SpherHarmonics::Ylm(0,0,theta,0.0)*hstar0;
	  Xlm10=ci*0.5*sqrt(2.0)*sqrt(12.0*PI)*SpherHarmonics::Ylm(1,0,theta,0.0)*hstar1;
	  Xlm11=ci*0.5*sqrt(2.0)*sqrt(12.0*PI)*SpherHarmonics::Ylm(1,1,theta,0.0)*hstar1;
      
	  // First do the case for S=0;
	  psi=psisymm;
	  // this refers to S=0, L=0, J=0 channel
	  psi+=Xlm00*(Misc::ceiphi(-2.0*delta_1s0)-1.0);
	  // S=0, L=0, J=0

	  psisquared=0.25*real(psi*conj(psi));
      
	  // Now let's do the case for S=1, M_S=+1
	  psi=psianti;
	  // S=1, L=1 and J=1,2,
	  psi+=Xlm10*(0.5*Misc::ceiphi(-2.0*delta_3p1)
		      +0.5*Misc::ceiphi(-2.0*delta_3p2)-1.0);
	  psia=Xlm11*0.5*(Misc::ceiphi(-2.0*delta_3p1)-Misc::ceiphi(-2.0*delta_3p2));
	  psisquared+=0.5*real(psi*conj(psi)+psia*conj(psia));
	  // Term is doubled to account for M_S=-1
      
	  // Now let's do the case with S=1, M_S=0;
	  psi=psianti;
	  // S=1, L=1, J=0,2
	  psi+=Xlm10*((2.0/3.0)*Misc::ceiphi(-2.0*delta_3p2)
		      +(1.0/3.0)*Misc::ceiphi(-2.0*delta_3p0)-1.0);
	  psia=Xlm11*(Misc::ceiphi(-2.0*delta_3p2)-Misc::ceiphi(-2.0*delta_3p0))/3.0;
	  psisquared+=0.25*real(psi*conj(psi)+2.0*psia*conj(psia));
	}

    }
  else psisquared=0.25*real(psisymm*conj(psisymm))
	 +0.75*real(psianti*conj(psianti));
  if(r>epsilon && psisquared<-0.001){
    printf("in CalcPsiSquared psisquared=%g,q=%g,r=%g\n",psisquared,q,r);
    if(r>epsilon) exit(1);
  }

  return	psisquared;
}

bool CWaveFunction_pp::GetPsiNoneCoupledState(long iq,
											   short n_L,
											   const char * szPartialWaveName,
											   short nTagPartialWave,
											   bool	bReidInput)
{
  const	double	ALPHA = 137.036;

	strcpy(pname, szPartialWaveName);

	// Attention!
	// Here the consistency of nTagPartialWave and szPartialWaveName is not checked, though they should match each other.

	double	fq = (iq + 0.5) * delq;
	double	fDelR;
	short	nq1q2 = 1;
	complex<double>	cPsi00Out, cPsi10Out, cPsi00In, cPsi10In, cPsi00, cPsi10, cPsi0, cPsi1, cPsi2;
	const	double	fm1 = MPROTON;
	const	double	fm2 = MPROTON;
	const	double  fmu= fm1 * fm2 / (fm1 + fm2);
	complex<double>	CI(0.0, 1.0);
	double	fx0, fx1;
	double	fE, fV;
	double	fEta;
	double	fR0, fR1, fR2;
	double	fSigma;

	fDelR = fRMax / double(nRMax);
	fE = fq * fq / (2.0 * fmu);
	fEta = fmu / (ALPHA * fq);

	fR1 = fRMax + 2.0 * fDelR;
	fR0 = fRMax + fDelR;

	fx1 = fR1 * fq / HBARC;
	fx0 = fR0 * fq / HBARC;

	//	Incoming
	using namespace	CoulWave;
	cPsi10In = CWincoming(n_L, fx1, fEta);
	cPsi00In = CWincoming(n_L, fx0, fEta);

	//	Outgoing
	cPsi10Out = conj(cPsi10In);
	cPsi00Out = conj(cPsi00In);

	// Get the phase to fit the boundary condition at r=0 point

	cPsi10 = cPsi10In;
	cPsi1 = cPsi10;
	cPsi00 = cPsi00In;
	cPsi0 = cPsi00;

	double	v11, v12, v22;

	if(true == bReidInput)
	{
		for(long i=nRMax;i>=0;i--)
		{
			fR2 = fR1;
			cPsi2 = cPsi1;
			fR1 = fR0;
			cPsi1 = cPsi0;
			fR0 = i * fDelR;
			fx0 = fR0 * fq / HBARC;

			double	fmr = sqrt(fmu * fmu + fq * fq);
			creid93(&fR1, pname, ptype, &v11, &v12, &v22);
			fV = nq1q2 * (HBARC / ALPHA)/ fR1  + n_L * (n_L + 1.0) * HBARC * HBARC / fR1 / fR1 / 2.0 / fmr;
			//Attention: In some books, they use the convention that U = -V and thus there is a sign change before V11, V22 and V12.
			fV += v11;
			cPsi0 = 2.0 * cPsi1 - cPsi2 + (2.0 * fmr * fV - fq * fq) * fDelR * fDelR * cPsi1 / (HBARC * HBARC);

			pcPsi[i] = cPsi0;
		}
	}
	else
	{
		for(long i=nRMax;i>=0;i--)
		{
			fR2 = fR1;
			cPsi2 = cPsi1;
			fR1 = fR0;
			cPsi1 = cPsi0;
			fR0 = i * fDelR;
			fx0 = fR0 * fq / HBARC;

			double	fmr = sqrt(fmu * fmu + fq * fq);
			fV = nq1q2 * (HBARC / ALPHA)/ fR1  + n_L * (n_L + 1.0)  * HBARC * HBARC / fR1 / fR1 / 2.0 / fmr;
			cPsi0 = 2.0 * cPsi1 - cPsi2 + (2.0 * fmr * fV - fq * fq) * fDelR * fDelR * cPsi1 / (HBARC * HBARC);

			pcPsi[i] = cPsi0;
		}
	}
	//	fPhase = atan2(imag(pcPsi[0]), real(pcPsi[0]));
	double	fPhase = 2.0 * atan2(imag(cPsi1), real(cPsi1)) - atan2(imag(cPsi2), real(cPsi2));

	//	Calculate the wave-function with the correct phase
	fSigma = fPhase;

	cPsi10 = (cPsi10In - cPsi10Out * exp(2.0 * CI * fSigma))/ 2.0;
	cPsi1 = cPsi10;
	cPsi00 = (cPsi00In - cPsi00Out * exp(2.0 * CI * fSigma))/ 2.0;
	cPsi0 = cPsi00;

	fR1 = fRMax + 2.0 * fDelR;
	fR0 = fRMax + fDelR;

	fx1 = fR1 * fq / HBARC;
	fx0 = fR0 * fq / HBARC;
	if(true == bReidInput)
	{
		for(long i=nRMax;i>=0;i--)
		{
			fR2 = fR1;
			cPsi2 = cPsi1;
			fR1 = fR0;
			cPsi1 = cPsi0;
			fR0 = i * fDelR;
			fx0 = fR0 * fq / HBARC;

			double	fmr = sqrt(fmu * fmu + fq * fq);
			creid93(&fR1, pname, ptype, &v11, &v12, &v22);
			fV = nq1q2 * (HBARC / ALPHA)/ fR1  + n_L * (n_L + 1.0)  * HBARC * HBARC / fR1 / fR1 / 2.0 / fmr;
			fV += v11;
			cPsi0 = 2.0 * cPsi1 - cPsi2 + (2.0 * fmr * fV - fq * fq) * fDelR * fDelR * cPsi1 / (HBARC * HBARC);

			pcPsi[i] = cPsi0;
		}
	}
	else
	{
		for(long i=nRMax;i>=0;i--)
		{
			fR2 = fR1;
			cPsi2 = cPsi1;
			fR1 = fR0;
			cPsi1 = cPsi0;
			fR0 = i * fDelR;
			fx0 = fR0 * fq / HBARC;

			double	fmr = sqrt(fmu * fmu + fq * fq);
			fV = nq1q2 * (HBARC / ALPHA)/ fR1  + n_L * (n_L + 1.0)  * HBARC * HBARC / fR1 / fR1 / 2.0 / fmr;
			cPsi0 = 2.0 * cPsi1 - cPsi2 + (2.0 * fmr * fV - fq * fq) * fDelR * fDelR * cPsi1 / (HBARC * HBARC);

			pcPsi[i] = cPsi0;
		}
	}

	switch(nTagPartialWave)
	{
	case	TAG_1S0:
		for(long ir=nRMax;ir>=0;ir--)
		{
			pc1S0[ir] = pcPsi[ir];
		}
		pfPhaseShift[0][iq] = fPhase;
		break;
	case	TAG_3P0:
		for(long ir=nRMax;ir>=0;ir--)
		{
			pc3P0[ir] = pcPsi[ir];
		}
		pfPhaseShift[1][iq] = fPhase;
		break;
	case	TAG_3P1:
		for(long ir=nRMax;ir>=0;ir--)
		{
			pc3P1[ir] = pcPsi[ir];
		}
		pfPhaseShift[2][iq] = fPhase;
		break;
	case	TAG_NS:
		for(long ir=nRMax;ir>=0;ir--)
		{
			pcNS[ir] = pcPsi[ir];
		}
		pfPhaseShiftNS[iq] = fPhase;
		break;
	case	TAG_3P2NOF:
		for(long ir=nRMax;ir>=0;ir--)
		{
			pc3P2NoF[ir] = pcPsi[ir];
		}
		pfPhaseShift3P2NoF[iq] = fPhase;
		break;
	default:
		break;
	}
	return true;
}
