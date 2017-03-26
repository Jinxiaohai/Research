#ifndef __INCLUDE_MISC_CC
#define __INCLUDE_MISC_CC
#include "misc.h"

bool Misc::comparestrings(char *s1,char *s2){
  int length1,length2,ic;
  bool answer;
  answer=0;
  length1=strlen(s1);
  length2=strlen(s2);
  if(length1==length2){
    answer=1;
    for(ic=0;ic<length1;ic++){
      if(s1[ic]!=s2[ic]){
	answer=0;
	goto NOMATCH;
      }
    }
  }
 NOMATCH:
  return answer;
}

double Misc::triangle(double m0,double m1,double m2){
  double answer,m0sq,m1sq,m2sq;
  if(m0<m1+m2) {
    printf("Disaster with triangle\n");
    exit(1);
  }
  m0sq=m0*m0;m1sq=m1*m1;m2sq=m2*m2;
  answer=m0sq*m0sq+m1sq*m1sq+m2sq*m2sq;
  answer=answer-2.0*(m0sq*m1sq+m0sq*m2sq+m1sq*m2sq);
  answer=answer/(4.0*m0sq);
  return answer;
}

double Misc::GetRapidity(double *pa){
  //double answer=0.5*log((pa[0]+pa[3])/(pa[0]-pa[3]));
  //if(answer>1.5) printf("answer is toooo big\n",answer);
  //return answer;
  return 0.5*log((pa[0]+pa[3])/(pa[0]-pa[3]));
}

double Misc::GetDely(double *pa,double *pb){
  return GetRapidity(pa)-GetRapidity(pb);
}

double Misc::GetQinv(double *pa,double *pb){
  int alpha;
  double answer;
  answer=0;
  answer=pow(pa[1]-pb[1],2)+pow(pa[2]-pb[2],2)+pow(pa[3]-pb[3],2)
    -pow(pa[0]-pb[0],2);
  if(answer<0.0){
    printf("DISASTER with GetQinv, wrong sign, = %g\n",answer);
    exit(1);
  }
  return sqrt(answer);
}

complex<double> Misc::cexp(complex<double> z){
  return exp(real(z))*ceiphi(imag(z));
}
complex<double> Misc::ceiphi(double phi){
  return complex<double>(cos(phi),sin(phi));
}

complex<double> Misc::cpow(complex<double> z,complex<double> a){
  complex<double> ci(0.0,1.0);
  double zr=real(z);
  double zi=imag(z);
  double phi=atan2(zi,zr);
  double zmag=sqrt(zr*zr+zi*zi);
  complex<double> alnz=a*(log(zmag)+ci*phi);
  return cexp(alnz);
}

void Misc::lorentz(double *u,double *p,double *pprime){
  const int n[4]={1,0,0,0};
  const int g[4]={1,-1,-1,-1};
  int mu;
  double udotn=0.0,pdotn=0.0,pdotu=0.0;
  //double mtest=0.0,mtestprime=0.0;
  for(mu=0;mu<4;mu++){
    //mtest=mtest+p[mu]*p[mu]*g[mu];
    pdotn=pdotn+p[mu]*n[mu]*g[mu];
    pdotu=pdotu+p[mu]*u[mu]*g[mu];
    udotn=udotn+u[mu]*n[mu]*g[mu];
  }
  for(mu=0;mu<4;mu++){
    pprime[mu]=-((pdotu+pdotn)/(1.0+udotn))*(n[mu]+u[mu])
      +2*pdotn*u[mu]+p[mu];
  }
}

int Misc::iround(double x){
  return int(floor(x+0.5));
}
// Using GSL
double Misc::cgc(double j1,double m1,double j2,double m2,double J,double M){
  const double e=1.0E-10;
  int sign=-1;
  if((lrint((M+j1-j2)))%2==0) sign=1;
  return sign*sqrt(2.0*J+1)
    *gsl_sf_coupling_3j(lrint(2*j1),lrint(2*j2),lrint(2*J),
			lrint(2*m1),lrint(2*m2),-lrint(2*M));
}

// from Alexander Volya, Declan Mulhall
//page 44 edmonds
double Misc::oldcgc(double j1,double m1,double j2,double m2,double j,double m){
  //rule check
  if (fabs(m1+m2-m)>0.01) return 0.0; //sum of m's must be zero
  if (j1+j2<j) return 0.0; //triangle rule
  if (j+j1<j2) return 0.0; 
  double p1,p2,p3,p4,p5,thesum = 0,ans;
  double z=0;
  if(m1 + m2 == m && m<=j && -j <=m){
    p1=(2*j+1)*cgc_fractorial((j1+j2-j),(j1+j2+j+1));
    p2=cgc_fractorial((j1-m1),(j1-j2+j));
    p3=cgc_fractorial((j2-m2),(j2-j1+j));
    p4=cgc_fractorial((j+m),(j1+m1));
    p5=cgc_fractorial((j-m),(j2+m2));
    for(z=0;z<=2*(j1+j2)+j;z++){
      if(j1-m1-z >=0 && j-m-z >=0 && j2-j+m1+z >=0 ){
	thesum = thesum +
	  pow(-1,z+j1-m1)*cgc_fractorial(j1+m1+z,j1-m1-z)
	  *cgc_fractorial(j2+j-m1-z,j-m-z)
	  /(cgc_factorial(z)*cgc_factorial(j2-j+m1+z));
      }
    }
    ans = sqrt(p1*p2*p3*p4*p5) * thesum;
    return ans;
  }
  else return 0;
}


double Misc::cgc_factorial(double n){
  double i,j=1;
  for (i=1;i<=n;i++) j=j*i;
  return j;
}

double Misc::cgc_fractorial(double n,double m){
  double  i,j=1;
  if(n>m){
    for(i=m+1;i<n+1;i++){j=j*i;}
    return j;
  }
  if(n<m){
    for(i=n+1;i<m+1;i++){j=j*i;}
    return 1/j;
  }
  return 1;
}

int Misc::cgc_delta (int x, int y) {if (x==y) return 1; else return 0;}

void Misc::outsidelong(double *pa,double *pb,
		 double &qinv,double &qout,double &qside,double &qlong){
  double vs,gamma,ptot[4],p[4],q[4],qtemp,ptot_perp;
  int alpha;
  for(alpha=0;alpha<4;alpha++){
    q[alpha]=0.5*(pa[alpha]-pb[alpha]);
    ptot[alpha]=pa[alpha]+pb[alpha];
  }

  // Perform long. comoving boost
  vs=ptot[3]/ptot[0];
  gamma=1.0/sqrt(1.0-vs*vs);
  ptot[0]=gamma*(ptot[0]-vs*ptot[3]);
  ptot[3]=0.0;
  qtemp=q[3];
  q[3]=gamma*(q[3]-vs*q[0]);
  q[0]=gamma*(q[0]-vs*qtemp);

  ptot_perp=sqrt(ptot[1]*ptot[1]+ptot[2]*ptot[2]);
  qout=(q[1]*ptot[1]+q[2]*ptot[2])/ptot_perp;
  qside=sqrt(q[2]*ptot[1]-q[1]*ptot[2])/ptot_perp;
  qlong=q[3];

  vs=ptot_perp/ptot[0];
  gamma=1.0/sqrt(1.0-vs*vs);
  qout=gamma*(qout-vs*q[0]);

}

#endif
