#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iostream>

using namespace std;

#include "coral.h"

int main(){
  int ir;
  double x,y,z,ctheta,phi,guess,r;
  const double PI=4.0*atan(1.0);
  CRandom randy(-12345);
  double Rx=3,Ry=4,Rz=5,Xoff=0.5,Yoff=1.5,Zoff=2.5;
  char sdirname[160];
  char ArrayParsFilename[120];
  CCHArray *Asource;
  CSourceCalc_Gaussian *scalc;

  // Create Array to store Cart. Harmonic source info
  sprintf(ArrayParsFilename,"parameters/apars_gauss.dat");
  Asource=new CCHArray(ArrayParsFilename);
  Asource->PrintPars();

  // Create Source Calc Object
  scalc=new CSourceCalc_Gaussian();
  scalc->SetSPars(Rx,Ry,Rz,Xoff,Yoff,Zoff);
  parameter::PrintPars(scalc->spars);

  // Have Source Calc Object Fill Array
  scalc->CalcS(Asource);
  scalc->NormCheck(Asource);
  scalc->CalcEffGaussPars(Asource);

  // Make 3D Cartesian array from Cart. Harmonic Data
  sprintf(ArrayParsFilename,"parameters/apars3d.dat\0");
  C3DArray *threed=new C3DArray(ArrayParsFilename);
  ArrayCalc::Calc3DArrayFromAExpArray(Asource,threed);
  printf("Created 3darray\n");

  // Calc S for some random angles and compare to correct answer
  for(ir=1;ir<Asource->GetNRADIAL();ir+=2){
    ctheta=1.0-2.0*randy.ran();
    phi=2.0*PI*randy.ran();
    r=(0.5+ir)*Asource->GetRADSTEP();
    z=r*ctheta;
    x=r*sqrt(1.0-ctheta*ctheta);
    y=x*sin(phi);
    x=x*cos(phi);
    guess=exp(-0.25*( ((x-Xoff)*(x-Xoff)/(Rx*Rx))
		      +((y-Yoff)*(y-Yoff)/(Ry*Ry))
		      +((z-Zoff)*(z-Zoff)/(Rz*Rz))));
    guess=guess/(Rx*Ry*Rz*pow(4.0*PI,1.5));
    printf("r=%g, S(%g,%g,%g)=%g =? %g =? %g\n",r,x,y,z,
	   Asource->AExpand(x,y,z),guess,
	   threed->GetElement_Interpolate(x,y,z));
  }
  
  // Write array info to file
  sprintf(sdirname,"sdata/gauss_R%g,%g,%g_off%g,%g,%g",
	  Rx,Ry,Rz,Xoff,Yoff,Zoff);
  Asource->WriteAllA(sdirname);
}

