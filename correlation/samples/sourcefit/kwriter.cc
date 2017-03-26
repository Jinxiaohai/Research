#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iostream>

using namespace std;

#include "coral.h"

int main(){
  CWaveFunction *wf;
  CKernel *kernel;
  char pardirname[80];
  sprintf(pardirname,"parameters/pp\0");

  char wfparsfilename[120];
  sprintf(wfparsfilename,"%s/kparameters.dat\0",pardirname);
  char kdatadirname[200];
  sprintf(kdatadirname,"kdata/pp\0");

  kernel=new CKernel(wfparsfilename);

  // Either read or calc
  //kernel->Read(kdatadirname);
  wf=new CWaveFunction_pp(wfparsfilename);
  kernel->Calc(wf);
  //kernel->Calc_ClassCoul(938.28,493.677,1);
  kernel->Write(kdatadirname);

}

