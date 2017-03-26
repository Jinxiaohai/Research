#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "random.cc"

int main(){
  CRandom random(-time(NULL));
  int i,nmc=1000000;
  double mass,pxbar,pybar,pzbar,pbar,ztest,T=100.0;
  double p[4];
  printf("Enter nmc : ");
  scanf("%d",&nmc);
  mass=1.0E-8;
  pxbar=pybar=pzbar=pbar=0.0;
  ztest=0.0;
  for(i=0;i<nmc;i++){
    random.generate_boltzmann(mass,T,p);
    pbar+=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
    pxbar+=p[1];
    pybar+=p[2];
    pzbar+=p[3];
    ztest+=p[1]*p[3];
  }
  pbar=pbar/double(nmc);
  ztest=ztest/double(nmc);
  printf("ztest=%g\n",ztest);
  printf("pbar=%g =? %g\n",pbar,3.0*T);
  pxbar=pxbar/double(nmc);
  pybar=pybar/double(nmc);
  pzbar=pzbar/double(nmc);
  printf("pxbar,pybar,pzbar=(%g,%g,%g) =? 0\n",pxbar,pybar,pzbar);
  return 0;
}
