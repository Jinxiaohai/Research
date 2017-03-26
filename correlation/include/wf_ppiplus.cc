#ifndef __CWAVEFUNCTION_WF_PN_CC
#define __CWAVEFUNCTION_WF_PN_CC

CWaveFunction_ppiplus::CWaveFunction_ppiplus(char *parsfilename) : CWaveFunction() {
  int iq,ichannel;
  double q;
  ParsInit(parsfilename);

  m1=MPROTON;
  m2=MPI;
  q1q2=1;
  if(COULOMB==0) q1q2=0;
  nchannels=3;
  ellmax=1;
  InitArrays();
  printf("Arrays Initialized\n");

  ell[0]=0;
  ell[1]=1;
  ell[2]=1;

  InitWaves();
  printf("Partial Waves Initialized\n");

  channelweight[0]=1.0;
  channelweight[1]=1.0; // J=1/2
  channelweight[2]=2.0; // J=3/2

  read_phaseshifts();

  for(ichannel=0;ichannel<nchannels;ichannel++){
    for(iq=0;iq<nqmax;iq++){
      q=qarray[iq];
      //printf("ichannel=%d, q=%g, delta=%g, ddeltadq=%g\n",
      //     ichannel,q,delta[ichannel][iq]*180.0/PI,
      //     ddeltadq[ichannel][iq]*180.0/PI);
      Wepsilon[ichannel][iq]=ddeltadq[ichannel][iq]
	-GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],delta[ichannel][iq])
	+GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],0.0);
      Wepsilon[ichannel][iq]=3.0*Wepsilon[ichannel][iq]
	/(4.0*PI*pow(epsilon,3));
    }
  }
  printf("Initialization finished\n");
}

double CWaveFunction_ppiplus::CalcPsiSquared(int iq,double r,double ctheta){
  double psisquared,x,dpsi2,q;
  double delta_s31,delta_p31,delta_p33;
  complex<double> psi,hstar,psi0;
  int ichannel;
  
  q=qarray[iq];
  if(iq>=nqmax){
    printf("iq too large!\n");
    exit(1);
  }
  psi0=planewave[iq]->planewave(r,ctheta);

  if(STRONG==1){
    if(r<epsilon){
      psisquared=real(psi0*conj(psi0));
      for(ichannel=0;ichannel<nchannels;ichannel++){
	dpsi2=channelweight[ichannel]*2.0*PI*Wepsilon[ichannel][iq]
	  *pow(HBARC,3)/(q*q);
	psisquared+=dpsi2;
      }
    }
    else{
      x=q*r/HBARC;
      delta_s31=delta[0][iq];
      delta_p31=delta[1][iq];
      delta_p33=delta[2][iq];
 
      psi=psi0;
      // this refers to m_s=+1/2 s wave
      ichannel=0;
      hstar=partwave[ell[ichannel]][iq]->GetPhiIncoming(r)/x;
      psi+=0.5*hstar*(Misc::ceiphi(-2.0*delta_s31)-1.0);
       
      // m_s still equals 1/2, but now ell=1
      ichannel=1;
      hstar=partwave[ell[ichannel]][iq]->GetPhiIncoming(r)/x;
      psi+=0.5*hstar*((2.0/3.0)*Misc::ceiphi(-2.0*delta_p33)
      		   +(1.0/3.0)*Misc::ceiphi(-2.0*delta_p31)-1.0)*ci*(3.0)*ctheta; 
      psisquared=real(psi*conj(psi));
       
      //Now consider m_s=-1/2, m_ell=1, Note: doesn't interfere with m_s=+1/2
      ichannel=2;
      psi=0.5*hstar*(Misc::ceiphi(-2.0*delta_p33)-Misc::ceiphi(-2.0*delta_p31))
	*ci*(3.0)*sqrt(1.0-ctheta*ctheta)/3.0;
      psisquared+=real(psi*conj(psi));

    }
  }
  else psisquared=real(psi0*conj(psi0));
  return psisquared;

}

void CWaveFunction_ppiplus::read_phaseshifts(){
#include "ppi_phaseshiftdat.cc"
  int iqdata,iq,ichannel;
  double w1,w2,q;
  for(ichannel=0;ichannel<3;ichannel++){
    for(iq=0;iq<nqmax;iq++){
      q=qarray[iq];
      iqdata=int(floor(q)/delqdata);
      if(iqdata<=nqdata){
	w1=(delqdata*double(iqdata+1)-q)/delqdata;
	w2=1.0-w1;
	delta[ichannel][iq]=w1*data_delta[ichannel][iqdata]
	  +w2*data_delta[ichannel][iqdata+1];
	ddeltadq[ichannel][iq]=w1*data_ddeltadq[ichannel][iqdata]
	  +w2*data_ddeltadq[ichannel][iqdata+1];
      }
      else{
	delta[ichannel][iq]=0.0;
	ddeltadq[ichannel][iq]=0.0;
      }
    }
  }

}

#endif
