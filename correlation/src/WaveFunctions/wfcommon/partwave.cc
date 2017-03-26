CPartWave::CPartWave(double etaset,int q1q2set,double qset,int ellset,
		     double epsilonset){
  epsilon=epsilonset;
  delr=0.01;
  nrmax=3000;
  complex<double> cg;
  PI=4.0*atan(1.0);
  ci=complex<double>(0.0,1.0);
  q=qset;
  q1q2=q1q2set;
  eta=etaset;
  ell=ellset;
  sigma=0.0;
  if(q1q2!=0){
    cg=CoulWave::cgamma(ell+1.0+ci*eta);
    sigma=atan2(imag(cg),real(cg));
  }
  phi_init();
}

CPartWave::CPartWave(double etaset,int q1q2set,double qset,int ellset,
		     double epsilonset,int nrmaxset,double delrset){
  epsilon=epsilonset;
  nrmax=nrmaxset;
  complex<double> cg;
  PI=4.0*atan(1.0);
  ci=complex<double>(0.0,1.0);
  q=qset;
  q1q2=q1q2set;
  eta=etaset;
  ell=ellset;
  sigma=0.0;
  if(q1q2!=0){
    cg=CoulWave::cgamma(ell+1.0+ci*eta);
    sigma=atan2(imag(cg),real(cg));
  }
  phi_init();
}

CPartWave::~CPartWave(){
  delete [] phi;
};

void CPartWave::phi_init(){
  const double HBARC=197.3269602;
  complex<double> phi0,phi1,phi2,phitest;
  int i,ir,ir0,ir00,oversample;
  double deltax,x,r;

  phi=new complex<double> [nrmax+1];

  if(q1q2!=0){
    for(ir=0;ir<nrmax;ir++){
      r=(0.5+ir)*delr;
      x=q*r/HBARC;
      phi[ir]=CoulWave::CWincoming(ell,x,eta)*Misc::ceiphi(-sigma);
    }
    
  }
  else{
    for(ir=0;ir<nrmax;ir++){
      r=(0.5+ir)*delr;
      x=q*r/HBARC;
      phi[ir]=x*Bessel::hstarn(ell,x);
    }
  }

}

complex<double> CPartWave::GetPhiIncoming(double r){
  const double HBARC=197.3269602;
  double x,a;
  complex<double> answer;
  int ir;
  ir=int(floor(r/delr));
  if(ir<nrmax){
    a=(r-delr*(ir+0.5))/delr;
    if(a>0.0 || ir==0){
      answer=(1.0-a)*phi[ir]+a*phi[ir+1];
    }
    else{
      answer=(1.0+a)*phi[ir]-a*phi[ir-1];
    }
  }
  else{
    x=q*r/HBARC;
    if(q1q2!=0){
      answer=CoulWave::CWincoming(ell,x,eta);
      answer=answer*Misc::ceiphi(-sigma);
    }
    else
      answer=x*Bessel::hstarn(ell,x);
  }
  return answer;
}

complex<double> CPartWave::GetPhiOutgoing(double r){
  complex<double> answer;
  answer=conj(GetPhiIncoming(r));
  if(q1q2!=0) answer=answer*Misc::ceiphi(-2.0*sigma);
  return answer;
}
