#ifndef __CWAVEFUNCTION_WF_PKPLUS_CC
#define __CWAVEFUNCTION_WF_PKPLUS_CC

CWaveFunction_pkplus::CWaveFunction_pkplus(char *parsfilename) : CWaveFunction(){
  int iq,ichannel;
  double q;
  ParsInit(parsfilename);

  m1=MPROTON;
  m2=MKAON;
  q1q2=1;
  if(COULOMB==0) q1q2=0;
  nchannels=3;
  ellmax=1;
  InitArrays();

  ell[0]=0;
  ell[1]=1;
  ell[2]=1;

  InitWaves();

  channelweight[0]=1.0;
  channelweight[1]=1.0;
  channelweight[2]=2.0;

  read_phaseshifts();

  for(ichannel=0;ichannel<nchannels;ichannel++){
    for(iq=0;iq<nqmax;iq++){
      //if(ichannel>=0){
      //delta[ichannel][iq]=ddeltadq[ichannel][iq]=0.0;
      //}
      q=qarray[iq];
      /* printf("ichannel=%d, q=%g, delta=%g, ddeltadq=%g\n",
	 ichannel,q,delta[ichannel][iq]*180.0/PI,
	 ddeltadq[ichannel][iq]*180.0/PI); */
      Wepsilon[ichannel][iq]=ddeltadq[ichannel][iq]
	-GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],delta[ichannel][iq])
	+GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],0.0);
      Wepsilon[ichannel][iq]=3.0*Wepsilon[ichannel][iq]
	/(4.0*PI*pow(epsilon,3));
    }
  }
}

double CWaveFunction_pkplus::CalcPsiSquared(int iq,double r,double ctheta){
  double q,psisquared,x,dpsi2;
  double delta_s11,delta_p11,delta_p13;
  complex<double> psi,hstar0,hstar1,psi0;
  int ipartial,ichannel;

  q=qarray[iq];
  x=q*r/HBARC;
  psi0=planewave[iq]->planewave(r,ctheta);
  if(iq>=nqmax){
    printf("iq too large!\n");
    exit(1);
  }

  if(r<epsilon){
    psisquared=real(psi0*conj(psi0));
    if(STRONG==1){
      for(ichannel=0;ichannel<nchannels;ichannel++){
	dpsi2=channelweight[ichannel]*2.0*PI*Wepsilon[ichannel][iq]
	  *pow(HBARC,3)/(q*q);
	psisquared+=dpsi2;
      }
    }
  }
  else{
    if(STRONG==1){
      hstar0=partwave[0][iq]->GetPhiIncoming(r)/x;
      hstar1=partwave[1][iq]->GetPhiIncoming(r)/x;

      delta_s11=delta[0][iq];
      delta_p11=delta[1][iq];
      delta_p13=delta[2][iq];
      
      // s-wave m_s=+1/2
      psi=psi0+0.5*hstar0*(Misc::ceiphi(-2.0*delta_s11)-1.0);
      
      // m_s still equals 1/2, but now ell=1
      psi+=0.5*hstar1*((2.0/3.0)*Misc::ceiphi(-2.0*delta_p13)
		       +(1.0/3.0)*Misc::ceiphi(-2.0*delta_p11)-1.0)*ci*(3.0)*ctheta;
      
      psisquared=real(psi*conj(psi));
      
      //Now consider m_s=-1/2, m_ell=1, Note: doesn't interfere with m_s=+1/2
      psi=0.5*hstar1*(Misc::ceiphi(-2.0*delta_p13)-Misc::ceiphi(-2.0*delta_p11))
	*ci*(3.0)*sqrt(1.0-ctheta*ctheta)/3.0;
      psisquared+=real(psi*conj(psi));
    }
    else{
      psi=psi0;
      psisquared=real(psi0*conj(psi0));
    }
  }
  if(psisquared>500.0){
    printf("_____________________________________________\n");
    printf("in CalcPsiSquared, q=%g,r=%g, psisquared=huge=%g\n",q,r,psisquared);
    printf("hstar0=(%g,%g), hstar1=(%g,%g),psi=(%g,%g)\n",real(hstar0),imag(hstar0),
	   real(hstar1),imag(hstar1),real(psi),imag(psi));
  }
  return psisquared;

}

void CWaveFunction_pkplus::read_phaseshifts(){
  int iq,iqdata,ichannel;
  double q,w1,w2;
  char dummy[200];
  FILE *fptr;
#include "pk_phaseshiftdat.cc"

  for(ichannel=0;ichannel<3;ichannel++){
    delta[ichannel][0]=ddeltadq[ichannel][0]=0.0;
    if(ichannel==0) ddeltadq[0][0]=data_delta[0][1]-data_delta[0][0];
    for(iq=1;iq<nqmax;iq++){
      q=qarray[iq];
      iqdata=int(floor(q)/delqdata);
      if(iqdata<nqdata){
	w1=(delqdata*double(iqdata+1)-q)/delqdata;
	w2=1.0-w1;
	delta[ichannel][iq]=w1*data_delta[ichannel][iqdata]
	  +w2*data_delta[ichannel][iqdata];
	ddeltadq[ichannel][iq]=w1*data_ddeltadq[ichannel][iqdata]
	  +w2*data_ddeltadq[ichannel][iqdata];
      }
      else{
	delta[ichannel][iq]=ddeltadq[ichannel][iq]=0.0;
      }
    }
  }
}

#endif
