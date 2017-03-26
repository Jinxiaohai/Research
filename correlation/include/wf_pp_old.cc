#ifndef __CWAVEFUNCTION_WF_PP_CC
#define __CWAVEFUNCTION_WF_PP_CC

CWaveFunction_pp_old::CWaveFunction_pp_old(char *parsfilename) : CWaveFunction() {
  int iq,ichannel;
  double q;
  ParsInit(parsfilename);

  m1=MPROTON;
  m2=m1;
  q1q2=1;
  if(COULOMB==0) q1q2=0;
  nchannels=4;

  ellmax=1;
  InitArrays();
  printf("Arrays Initialized\n");
  ell[0]=0;
  ell[1]=ell[2]=ell[3]=1;

  InitWaves();
  printf("Partial Waves Initialized\n");

  // Channel weight is (2J+1)/[(2s1+1)*(2s2+1)]
  channelweight[0]=2*0.25;
  channelweight[1]=2*0.25;
  channelweight[2]=2*0.75;
  channelweight[3]=2*1.25;
  read_phaseshifts();
  printf("phaseshifts read in\n");

  for(ichannel=0;ichannel<nchannels;ichannel++){
    for(iq=0;iq<nqmax;iq++){
      q=qarray[iq];
      /*if(ichannel==0) printf("q=%g, ddeltadq=%g, W=%g, W0=%g\n",
	qarray[iq],ddeltadq[ichannel][iq],
	GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],delta[ichannel][iq]),
	GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],0.0));*/
      Wepsilon[ichannel][iq]=ddeltadq[ichannel][iq]
	-GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],delta[ichannel][iq])
	+GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],0.0);
      Wepsilon[ichannel][iq]=3.0*Wepsilon[ichannel][iq]
	*pow(HBARC/epsilon,3)/(2.0*q*q);
      /*if(ichannel==0){
	printf("q=%g, W=%g, wprime=%g\n",q,Wepsilon[ichannel][iq],
	ddeltadq[ichannel][iq]*3.0*pow(HBARC/epsilon,3)/(2.0*q*q));
	}
	else{
	Wepsilon[ichannel][iq]=0.0;
	delta[ichannel][iq]=0.0;
	}*/
    }
  }
  printf("Initialization finished\n");
}

double CWaveFunction_pp_old::CalcPsiSquared(int iq,double r,double ctheta){
  double psisquared,x,dpsi2,q,theta;
  double delta_1s0,delta_3s1,delta_3p0,delta_1p1,delta_3p1,delta_3p2;
  complex<double> psi,psisymm,psianti,psia,psib,hstar0,hstar1;
  complex<double> Xlm00,Xlm10,Xlm11;
  int ichannel;

  q=qarray[iq];
  if(iq>=nqmax){
    printf("iq too large!\n");
    exit(1);
  }
  psia=planewave[iq]->planewave(r,ctheta);
  psib=planewave[iq]->planewave(r,-ctheta);
  psisymm=(1.0/sqrt(2.0))*(psia+psib);
  psianti=(1.0/sqrt(2.0))*(psia-psib);
  
  if(STRONG==1){
    if(r<epsilon){
      psisquared=0.25*real(psisymm*conj(psisymm))
	+0.75*real(psianti*conj(psianti));
      for(ichannel=0;ichannel<nchannels;ichannel++){
	dpsi2=channelweight[ichannel]*Wepsilon[ichannel][iq];
	psisquared+=dpsi2;
      }
    }
    else{
      theta=acos(ctheta);
      x=q*r/HBARC;
      // Notation is (2S+1)-L-J
      delta_1s0=delta[0][iq];
      delta_3p0=delta[1][iq];
      delta_3p1=delta[2][iq];
      delta_3p2=delta[3][iq];
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
      // S=0, L=1, J=1
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
  //if(psisquared<0.0){
  //  printf("psisquared <0, = %g, r=%g, q=%g\n",psisquared,r,q);
  //}
  return psisquared;

}

void CWaveFunction_pp_old::read_phaseshifts(){
#include "pp_phaseshiftdat.cc"
  int iqdata,iq,ichannel;
  double w1,w2,q;
  for(ichannel=0;ichannel<4;ichannel++){
    for(iq=0;iq<nqmax;iq++){
      q=qarray[iq];
      iqdata=int(floor(q)/delqdata);
      if(iqdata<nqdata){
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
      //if(ichannel>0) delta[ichannel][iq]=ddeltadq[ichannel][iq]=0.0;
    }
  }
}

#endif
