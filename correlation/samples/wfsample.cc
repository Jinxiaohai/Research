#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iostream>

using namespace std;

#include "coral.h"

int main(){
  CWaveFunction *wf;
  double q,r,ctheta,R=3.0,x,y,z;
  int iq,Nq;
  double *c;
  int imc,NMC=100000;
  CRandom random(-256);
  char parsfilename[120];

  printf("Enter N_MonteCarlo : ");
  scanf("%d",&NMC);

  sprintf(parsfilename,"parameters/wfparameters.dat\0");
  wf=new CWaveFunction_pp(parsfilename);
  //wf->PrintCdelta(R,R,R);

  Nq=wf->GetNQMAX();
  c=new double[Nq];
  for(iq=0;iq<Nq;iq++){
    q=wf->GetQ(iq);
    c[iq]=0.0;
    for(imc=0;imc<NMC;imc++){
      x=R*sqrt(2.0)*random.gauss();
      y=R*sqrt(2.0)*random.gauss();
      z=R*sqrt(2.0)*random.gauss();
      r=sqrt(x*x+y*y+z*z);
      ctheta=z/r;
      //c[iq]+=wf->GetPsiSquared(q,r,ctheta);
      c[iq]+=wf->CalcPsiSquared(iq,r,ctheta);
    }
    c[iq]=c[iq]/double(NMC);
    printf("%5.2f : %g\n",q,c[iq]);
  }

}

