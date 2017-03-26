#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iostream>

using namespace std;

#include "coral.h"

int main(){
  CWaveFunction *wf;
  CKernel *kernel;
  double q,r,R=3.0,x,y,z,kvalue;
  int iq,Nq;
  double *c;
  int imc,NMC=100000;
  CRandom random(-123);
  char parsfilename[120];
  char kdatadirname[200];
  sprintf(kdatadirname,"kdata/pkplus");

  //printf("Enter N_MonteCarlo : ");
  //scanf("%d",&NMC);
  NMC=10000;

  sprintf(parsfilename,"parameters/wfparameters.dat\0");
  wf=new CWaveFunction_pkplus(parsfilename);
  kernel=new CKernel(parsfilename);
  // Either read or calc
  //kernel->Read(kdatadirname);
  kernel->Calc(wf);
  //kernel->Calc_ClassCoul(938.28,493.677,1);

  //kernel->Print();
  kernel->Write(kdatadirname);

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
      kvalue=kernel->GetValue(0,q,r);
      c[iq]+=1.0+kernel->GetValue(0,q,r);
      if(kvalue<-1.0 && r>1.0){
	printf("WARNING: r=%g, kernel BELOW ZERO, K=%g!!!!!!!!\n",r,kvalue);
	// exit(1);
      }
    }
    c[iq]=c[iq]/double(NMC);
    printf("c(q=%5.2f) = %g\n",q,c[iq]);
  }

  int ir;
  double ctheta;
  for(iq=0;iq<Nq;iq++){
    q=wf->GetQ(iq);
    for(ir=0;ir<kernel->GetNRMAX();ir++){
      ctheta=1.0-2.0*random.ran();
      r=(0.5+ir)*kernel->GetDELR();
      printf("iq=%d ir=%d ctheta=%g, ",iq,ir,ctheta);
      printf("phi^2=%g =? %g\n",
	     kernel->GetPsiSquared(q,r,ctheta),wf->CalcPsiSquared(iq,r,ctheta));
    }
  }
}

