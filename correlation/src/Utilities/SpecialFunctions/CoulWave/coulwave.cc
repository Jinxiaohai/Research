#ifndef __COULWAVE_CC
#define __COULWAVE_CC
#define NO_GSLCOULWAVE_BUG

using namespace std;

complex<double> CoulWave::cgamma(complex<double> z){
  complex<double> ci(0.0,1.0);
  gsl_sf_result gsl_lnr,gsl_arg;
  int dummy=gsl_sf_lngamma_complex_e(real(z),imag(z),&gsl_lnr,&gsl_arg);
  return exp(gsl_lnr.val)*(cos(gsl_arg.val)+ci*sin(gsl_arg.val));
}

/* ************************************* */

void CoulWave::GetFG(int L,double x,double eta,double *FL,double *GL){
  const double PI=4.0*atan(1.0);
  double expF,expG,phase,sigma;
  complex<double> cg,ci(0.0,1.0);
  double *fc,*gc;
  int k=0;
  fc=new double[L+1];
  gc=new double[L+1];
  cg=CoulWave::cgamma(double(L+1)+ci*eta);
  sigma=atan2(imag(cg),real(cg));

  // This calculates fc and gc arrays for indices L to L+k  
  int dummy=gsl_sf_coulomb_wave_FG_array(0,L,eta,x,fc,gc,&expF,&expG);
  *FL=fc[L]*exp(expF);
  *GL=gc[L]*exp(expG);
#ifndef NO_GSLCOULWAVE_BUG
  if(x>2.0){
    phase=atan2(*GL,*FL);
    phase=phase+x-0.5*PI*(L+1)+sigma-eta*log(2.0*x);
    phase+=0.5*double(L*(L+1))/x;
    phase=phase-2.0*PI*floor(0.5*phase/PI);
    if(phase>PI) phase=phase-2.0*PI;
    if(fabs(phase)>2.0){
      *FL=-*FL;
      *GL=-*GL;
      //printf("L=%d, x=%6.3f, phase=%10.6f, sigma=%g\n",L,x,phase,sigma);
    }
  }
#endif
  delete [] fc;
  delete [] gc;
}

// This calcules dCW/dx
void CoulWave::GetFGprime(int L,double x,double eta,double *FL,double *GL,
		      double *FLprime,double *GLprime){
  double expF,expG;
  double *fc,*gc,*fcp,*gcp;
  int k=0;
  fc=new double[k+1];
  gc=new double[k+1];
  fcp=new double[k+1];
  gcp=new double[k+1];
  // This calculates fc and gc arrays for indices L to L+k  
  int dummy=gsl_sf_coulomb_wave_FGp_array(L,k,eta,x,fc,fcp,gc,gcp,
					 &expF,&expG);
  *FL=fc[0]*exp(expF);
  *GL=gc[0]*exp(expG);
  *FLprime=fcp[0]*exp(expF);
  *GLprime=gcp[0]*exp(expG);
  delete [] fc;
  delete [] gc;
  delete [] fcp;
  delete [] gcp;
}

complex<double> CoulWave::CWoutgoing(int ell,double x,double eta){
  double FL,GL;
  complex<double> ci(0.0,1.0);
  GetFG(ell,x,eta,&FL,&GL);
  return FL-ci*GL;
} 

complex<double> CoulWave::CWincoming(int ell,double x,double eta){
  double FL,GL;
  complex<double> ci(0.0,1.0);
  GetFG(ell,x,eta,&FL,&GL);
  return FL+ci*GL;
}

/* ************************************* */

void CoulWave::phaseshift_CoulombCorrect(int ell,double q,double eta,
					 double &delta,double &ddeltadq){
  const double PI=4.0*atan(1.0);
  double *gamow,*dgamowdq;
  double tandelta0,tandelta,dtandelta0dq,dtandeltadq,x,y;
  int i;

  if(fabs(eta)>1.0E-10){
    gamow=new double[ell+1];
    dgamowdq=new double[ell+1];
    gamow[0]=2.0*PI*eta/(exp(2.0*PI*eta)-1.0);
    dgamowdq[0]=(gamow[0]/q)*(gamow[0]*exp(2.0*PI*eta)-1.0);
    for(i=0;i<ell;i++){
      gamow[i+1]=gamow[i]*(double((i+1)*(i+1))+eta*eta);
      dgamowdq[i+1]=dgamowdq[i]*(double((i+1)*(i+1))+eta*eta)
	-(2.0*eta*eta/q)*gamow[i];
    }
    
    tandelta0=tan(delta);
    dtandelta0dq=ddeltadq*(1.0+tandelta0*tandelta0);
    
    tandelta=tandelta0*gamow[ell];
    if(cos(delta)>0.0) delta=atan(tandelta);
    else delta=atan(tandelta)+PI;
    if(delta>PI) delta=delta-2.0*PI;
    
    dtandeltadq=gamow[ell]*dtandelta0dq+tandelta0*dgamowdq[ell];
    ddeltadq=dtandeltadq/(1.0+tandelta*tandelta);
    delete [] gamow;
    delete [] dgamowdq;
  }
}

#endif
