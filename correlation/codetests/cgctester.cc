#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <complex>

#include "misc.cc"
#include "random.cc"

using namespace std;
using namespace Misc;

int main(){
  CRandom randy(12345);
  double J,j1,j2,m1,m2,M;
  int i,n;
  double ans,oldans;
  
  printf("Enter NMC : ");
  scanf("%d",&n);
  //n=50;
  for(i=0;i<n;i++){
    j1=randy.iran(31);
    j2=randy.iran(31);
    j1*=0.5;    j2*=0.5;
    m1=j1-randy.iran(int(2*j1+1));
    m2=j2-randy.iran(int(2*j2+1));
    M=m1+m2;
    J=fabs(j1-j2);
    if(j1>j2){
      J+=randy.iran(int(j2+1));
    }
    else  J+=randy.iran(int(j1+1));
    if(J>(j1+j2+0.00001) || J<fabs(j1-j2)-0.00000001){
      printf("OUCH, J=%g, j1=%g, j2=%g\n",J,j1,j2);
      exit(1);
    }

    ans=cgc(j1,m1,j2,m2,J,M);
    oldans=oldcgc(j1,m1,j2,m2,J,M);
    if(fabs(ans-oldans)>1.0E-5){
      printf("J=%g, j1=%g, m1=%g, j2=%g, m2=%g,         %g =? %g\n",
	     J,j1,m1,j2,m2,ans,oldans);
      exit(1);
    }

    //printf("J=%g, j1=%g, m1=%g, j2=%g, m2=%g,         %g =? %g\n",
    //	   J,j1,m1,j2,m2,ans,oldans);

  }
  return 0;
  
}

