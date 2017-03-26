#include <cstdlib>
#include <cmath>
#include <cstdio>

using namespace std;

#include "coral.h"

int main(){
  int ir,nsample;
  char OSCARfilename[160],sdirname[160];
  char ArrayParsFilename[160];
  CCHArray *Asource;

  // Create Array for storing source info
  sprintf(ArrayParsFilename,"parameters/apars_oscar.dat");
  Asource=new CCHArray(ArrayParsFilename);

  // Initialize Source Calc Object
  CSourceCalc_OSCAR *scalc;
  scalc=new CSourceCalc_OSCAR();
  parameter::set(scalc->spars,"OSCARfilename",
		 "../../gromit_hbt/output_data/bjoscar_pionsonly.output.tmp");
  parameter::set(scalc->spars,"Pt",100);
  scalc->SetIDs(211,211);
  parameter::PrintPars(scalc->spars);

  // Calculate source array
  scalc->CalcS(Asource);
  scalc->NormCheck(Asource);
  scalc->CalcEffGaussPars(Asource);

  // Write source info to file
  sprintf(sdirname,"sdata/OSCAR_pionsonly_kt%g\0",
	  parameter::getD(scalc->spars,"Pt",-200)/2);
  Asource->WriteAllA(sdirname);
}

