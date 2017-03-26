#ifndef __INCLUDE_SOURCECALC_GX1D
#define __INCLUDE_SOURCECALC_GX1D

using namespace std;

CSourceCalc_GX1D::CSourceCalc_GX1D(){
  InitSPars();
}

void CSourceCalc_GX1D::InitSPars(){
  // DEFAULT VALUES
  parameter::set(spars,"lambdaG",0.3);
  parameter::set(spars,"lambdaX",0.3);
  parameter::set(spars,"R",5.0);
  parameter::set(spars,"X",10.0);
  parameter::set(spars,"a",5.0);
}

void CSourceCalc_GX1D::SetSPars(double lambdaGset,double Rset,
				double lambdaXset,double Xset,double aset){
  parameter::set(spars,"lambdaG",lambdaGset);
  parameter::set(spars,"R",Rset);
  parameter::set(spars,"lambdaX",lambdaXset);
  parameter::set(spars,"X",Xset);
  parameter::set(spars,"a",aset);
}

void CSourceCalc_GX1D::CalcS(int lx,int ly,int lz,CCHArray *A){
  const double PI=4.0*atan(1.0);

  int ir,L;
  int nradial;
  double lambdaG,lambdaX,R,X,a;
  double r,z,delr,normX,normG,value;
  lambdaG=parameter::getD(spars,"lambdaG",0.3);
  R=parameter::getD(spars,"R",5.0);
  lambdaX=parameter::getD(spars,"lambdaX",0.3);
  X=parameter::getD(spars,"X",10.0);
  a=parameter::getD(spars,"a",5.0);
  normG=pow(4.0*R*R*PI,-1.5);
  z=2.0*R*R/(X*X);
  normX=1.0/(4.0*PI*X*X*X*(2.0*z*Bessel::K1(z)+z*z*Bessel::K0(z)));

  nradial=A->GetNRADIAL();
  delr=A->GetRADSTEP();

  L=lx+ly+lz;

  for(ir=0;ir<nradial;ir++){
    r=delr*(0.5+ir);
    value=lambdaX*normX*exp(-sqrt(((r*r)/(X*X))+z*z));
    value+=lambdaG*normG*exp(-0.25*r*r/(R*R));
    if(L>0) value*=pow(r*r/(r*r+a*a),0.5*double(L));
    A->SetElement(lx,ly,lz,ir,value);
  }

}

void CSourceCalc_GX1D::CalcS(CCHArray *A){
  int dlx=1,dly=1,dlz=1;
  int lx,ly,lz,lmax;
  lmax=A->GetLMAX();
  if(A->GetXSYM()) dlx=2;
  if(A->GetYSYM()) dly=2;
  if(A->GetZSYM()) dlz=2;
  for(lx=0;lx<2;lx+=dlx){
    for(ly=0;ly<lmax-lx;ly+=dly){
      for(lz=0;lz<lmax-lx-ly;lz+=dlz){
	CalcS(lx,ly,lz,A);
      }
    }
  }
}

#endif
