#ifndef SCOTTUTILS_YLM_CC
#define SCOTTUTILS_YLM_CC

using namespace std;

double SpherHarmonics::legendre(int ell,double ctheta){
  return gsl_sf_legendre_Pl(ell,ctheta);
}

complex<double> SpherHarmonics::Ylm(int ell,int m,double theta,double phi){
  double ctheta;
  complex<double> answer;
  complex<double> ci(0.0,1.0);
  ctheta=cos(theta);
  answer=gsl_sf_legendre_sphPlm(ell,abs(m),ctheta)*Misc::ceiphi(m*phi);
  if(m<0) answer=answer*pow(-1.0,abs(m));
  return answer;
}
 
#endif
