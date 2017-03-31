#ifndef __INCLUDE_CHARRAY_CC
#define __INCLUDE_CHARRAY_CC

using namespace std;

CRandom *CCHArray::randy=NULL;
CCHCalc *CCHArray::chcalc=NULL;

CCHArray::CCHArray(int LMAXset,int NRADIALset,double RADSTEPset){
  if(chcalc==NULL) chcalc=new CCHCalc();
  NRADIAL=NRADIALset;
  LMAX=LMAXset;
  RADSTEP=RADSTEPset;
  XSYM=YSYM=ZSYM=0;
  dlx=dly=dlz=1;
  CreateArray();
}

CCHArray::CCHArray(int LMAXset,int NRADIALset,double RADSTEPset,
			       bool XSYMset,bool YSYMset,bool ZSYMset){
  if(chcalc==NULL) chcalc=new CCHCalc();
  NRADIAL=NRADIALset;
  RADSTEP=RADSTEPset;
  LMAX=LMAXset;
  XSYM=XSYMset;
  YSYM=YSYMset;
  ZSYM=ZSYMset;
  dlx=dly=dlz=1;
  if(XSYM) dlx=2;
  if(YSYM) dly=2;
  if(ZSYM) dlz=2;
  CreateArray();
}

CCHArray::CCHArray(char *arrayparsfilename){
  dlx=dly=dlz=1;
  if(chcalc==NULL) chcalc=new CCHCalc();
  parametermap apars;

  parameter::set(apars,"NRADIAL",25);
  parameter::set(apars,"RADSTEP",1.0);
  parameter::set(apars,"LMAX",4);
  // NRADIAL, IDENTICAL and RADSTEP are only used if wf is used by kernel 
  parameter::set(apars,"IDENTICAL",bool(0));
  parameter::set(apars,"XSYM",bool(0));
  parameter::set(apars,"YSYM",bool(0));
  parameter::set(apars,"ZSYM",bool(0));

  parameter::ReadParsFromFile(apars,arrayparsfilename);
  if(parameter::getB(apars,"IDENTICAL",0)){
    parameter::set(apars,"XSYM",bool(1));
    parameter::set(apars,"YSYM",bool(1));
    parameter::set(apars,"ZSYM",bool(1));
  }

  NRADIAL=parameter::getI(apars,"NRADIAL",-999);
  LMAX=parameter::getI(apars,"LMAX",-999);
  RADSTEP=parameter::getD(apars,"RADSTEP",-999);
  XSYM=parameter::getB(apars,"XSYM",0);
  YSYM=parameter::getB(apars,"YSYM",0);
  ZSYM=parameter::getB(apars,"ZSYM",0);
  if(XSYM) dlx=2;
  if(YSYM) dly=2;
  if(ZSYM) dlz=2;

  CreateArray();
}


void CCHArray::CreateArray(){
  int lx,ly,lz,ir;
  //PrintPars();
  A=new double ***[LMAX+1];
  for(lx=0;lx<=LMAX;lx+=dlx){
    A[lx]=new double **[LMAX+1-lx];
    for(ly=0;ly<=LMAX-lx;ly+=dly){
      A[lx][ly]=new double *[LMAX+1-lx-ly];
      for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
	A[lx][ly][lz]=new double[NRADIAL];
	for(ir=0;ir<NRADIAL;ir++) A[lx][ly][lz][ir]=0.0;
      }
    }
  }
  
}

CCHArray::~CCHArray(){
  int lx,ly,lz;
  for(lx=0;lx<=LMAX;lx+=dlx){
    for(ly=0;ly<=LMAX-lx;ly+=dly){
      for(lz=0;lz<=LMAX-lx-ly;lz+=dlz) delete [] A[lx][ly][lz];
      delete [] A[lx][ly];
    }
    delete [] A[lx];
  }
  delete A;
}

int CCHArray::GetLMAX(){
  return LMAX;
}

void CCHArray::SetLMAX(int LMAXset){
  if(LMAXset<LMAX) LMAX=LMAXset;
}

int CCHArray::GetNRADIAL(){
  return NRADIAL;
}

double CCHArray::GetRADSTEP(){
  return RADSTEP;
}

void CCHArray::SetRADSTEP(double RADSTEPset){
  RADSTEP=RADSTEPset;
}

double CCHArray::GetElement(int lx,int ly,int lz,int ir){
  int L=lx+ly+lz;
  if(L<=LMAX && lx%dlx==0 && ly%dly==0 && lz%dlz==0 && ir<NRADIAL) 
    return A[lx][ly][lz][ir];
  else{
    printf("WARNING: Tried to get non-existent element in CCHArray\n");
    printf("LMAX=%d, XSYM=%d, YSYM=%d, ZSYM=%d, ir=%d\n",
	   LMAX,XSYM,YSYM,ZSYM,ir);
    printf("l=(%d,%d,%d), ir=%d\n",lx,ly,lz,ir);
    return 0.0;
  }
}

double CCHArray::GetElement(int lx,int ly,int lz,double r){
  int ir=int(floor(r/RADSTEP));
  return GetElement(lx,ly,lz,ir);
}

void CCHArray::SetElement(int lx,int ly,int lz,int ir,double Element){
  int L=lx+ly+lz;
  if(L<=LMAX && L>=0 && lx%dlx==0 && ly%dly==0 && lz%dlz==0 && ir<NRADIAL) 
    A[lx][ly][lz][ir]=Element;
  else{
    printf("Tried to set element out of bounds in CCHArray\n");
    printf("LMAX=%d, XSYM=%d, YSYM=%d, ZSYM=%d, ir=%d\n",
	   LMAX,XSYM,YSYM,ZSYM,ir);
    printf("l=(%d,%d,%d), ir=%d\n",lx,ly,lz,ir);
    exit(1);
    
  }
}

void CCHArray::SetElement(int lx,int ly,int lz,double r,double Element){
  int ir=int(floor(r/RADSTEP));
  SetElement(lx,ly,lz,ir,Element);
}

void CCHArray::IncrementElement(int lx,int ly,int lz,int ir,double increment){
  int L=lx+ly+lz;
  if(L<=LMAX && L>=0 && lx%dlx==0 && ly%dly==0 && lz%dlz==0 && ir<NRADIAL) {
    A[lx][ly][lz][ir]+=increment;
  }
  else{
    printf("OUT OF BOUNDS IN INCREMENTELEMENT\n");
    printf("ell=(%d,%d,%d), ir=%d\n",lx,ly,lz,ir);
    exit(1);
  }

}

void CCHArray::IncrementElement(int lx,int ly,int lz,double r,double increment){
  int ir=int(floor(r/RADSTEP));
  IncrementElement(lx,ly,lz,ir,increment);
}

void CCHArray::ScaleArray(double scalefactor){
  int lx,ly,lz,ir;
    for(lx=0;lx<=LMAX;lx+=dlx){
      for(ly=0;ly<=LMAX-lx;ly+=dly){
	for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
	  for(ir=0;ir<NRADIAL;ir++) A[lx][ly][lz][ir]*=scalefactor;
      }
    }
  }
}

void CCHArray::ScaleArray(double scalefactor, int ir){
  int lx,ly,lz;
    for(lx=0;lx<=LMAX;lx+=dlx){
      for(ly=0;ly<=LMAX-lx;ly+=dly){
	for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
	  A[lx][ly][lz][ir]*=scalefactor;
      }
    }
  }
}

void CCHArray::ZeroArray_Partial(int LMAX_partial){
  int lx,ly,lz,ir;
    for(lx=0;lx<=LMAX;lx+=dlx){
      for(ly=0;ly<=LMAX-lx-ly;ly+=dly){
	for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
	  for(ir=0;ir<NRADIAL;ir++) A[lx][ly][lz][ir]=0.0;
      }
    }
  }
}

void CCHArray::ZeroArray_Partial(int LMAX_partial,int ir){
  int lx,ly,lz;
    for(lx=0;lx<=LMAX;lx+=dlx){
      for(ly=0;ly<=LMAX-lx-ly;ly+=dly){
	for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
	  A[lx][ly][lz][ir]=0.0;
      }
    }
  }
}

void CCHArray::ZeroArray(){
  int lx,ly,lz,ir;
  for(lx=0;lx<=LMAX;lx+=dlx){
    for(ly=0;ly<=LMAX-lx;ly+=dly){
      for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
	for(ir=0;ir<NRADIAL;ir++) A[lx][ly][lz][ir]=0.0;
      }
    }
  }
}

void CCHArray::ZeroArray(int lx,int ly,int lz){
  int ir;
  for(ir=0;ir<NRADIAL;ir++) A[lx][ly][lz][ir]=0.0;
}

void CCHArray::ZeroArray(int ir){
  int lx,ly,lz;
  for(lx=0;lx<=LMAX;lx+=dlx){
    for(ly=0;ly<=LMAX-lx;ly+=dly){
      for(lz=0;lz<=LMAX-lx-ly;lz+=dlz)
	A[lx][ly][lz][ir]=0.0;
    }
  }
}

void CCHArray::PrintArrayFixedIR(int ir){
  PrintArrayFixedIR(LMAX,ir);
}
void CCHArray::PrintArrayFixedIR(int LMAXPrint,int ir){
  int L,lx,ly,lz,dL=1;
  if(XSYM && YSYM && ZSYM) dL=2;
  if(LMAXPrint>LMAX) LMAXPrint=LMAX;
  
  printf("\n______________________________________________________\n",ir);
  for(L=0;L<=LMAXPrint;L+=dL){
    printf("     L=%d\n",L);
    printf(" lx\\ly:");
    for(ly=0;ly<=L;ly+=dly) printf(" %4d       ",ly);
    printf("\n");
      
    for(lx=0;lx<=L;lx+=dlx){
      printf(" %3d ",lx);
      for(ly=0;ly<=L-lx;ly+=dly){
	lz=L-lx-ly;
	if(ZSYM==0 || lz%2==0)
	  printf(" %10.3e ",A[lx][ly][lz][ir]);
	else printf("%10.3e ",0.0);
      }
      printf("\n");
    }
    printf("_________________________________________\n");
  }
}

double CCHArray::GetBiggest(int ir){
  int lx,ly,lz;
  double biggy=0.0;
  for(lx=0;lx<=LMAX;lx+=dlx){
    for(ly=0;ly<=LMAX-lx;ly+=dly){
      for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
	if(fabs(A[lx][ly][lz][ir])>fabs(biggy))
	  biggy=fabs(A[lx][ly][lz][ir]);
      }
    } 
  }
  return biggy;
}

bool CCHArray::GetXSYM(){
  if(XSYM) return 1;
  else return 0;
}

bool CCHArray::GetYSYM(){
  if(YSYM) return 1;
  else return 0;
}

bool CCHArray::GetZSYM(){
  if(ZSYM) return 1;
  else return 0;
}

void CCHArray::WriteAX(char *dirname){
  char filename[160],shellcommand[200];
  int ir,lx,ly,lz;
  FILE *fptr;

  sprintf(shellcommand,"mkdir -p %s\0",dirname);
  system(shellcommand);

  for(lx=0;lx<=1;lx+=dlx){
    for(ly=0;ly<=LMAX-lx;ly+=dly){
      for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
	sprintf(filename,"%s/lx%d_ly%d_lz%d.tmp\0",dirname,lx,ly,lz);
	fptr=fopen(filename,"w");
	fprintf(fptr,"%d %g\n",NRADIAL,RADSTEP);
	for(ir=0;ir<NRADIAL;ir++){
	  fprintf(fptr,"%g\n",A[lx][ly][lz][ir]);
	}
	fclose(fptr);
      }
    }
  }
}

void CCHArray::WriteAllA(char *dirname){
  char filename[160],shellcommand[200];
  int ir,lx,ly,lz;
  FILE *fptr;

  sprintf(shellcommand,"mkdir -p %s\0",dirname);
  system(shellcommand);

  for(lx=0;lx<=LMAX;lx+=dlx){
    for(ly=0;ly<=LMAX-lx;ly+=dly){
      for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
	sprintf(filename,"%s/lx%d_ly%d_lz%d.tmp\0",dirname,lx,ly,lz);
	fptr=fopen(filename,"w");
	fprintf(fptr,"%d %g\n",NRADIAL,RADSTEP);
	for(ir=0;ir<NRADIAL;ir++){
	  fprintf(fptr,"%16.10e\n",A[lx][ly][lz][ir]);
	}
	fclose(fptr);
      }
    }
  }
}

void CCHArray::ReadAX(char *dirname){
  char filename[160],dummy[200],shellcommand[200];
  int ir,lx,ly,lz,NRADIALread;
  double aa,RADSTEPread;
  FILE *fptr;

  for(lx=0;lx<=1;lx+=dlx){
    for(ly=0;ly<=LMAX-lx;ly+=dly){
      for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
	sprintf(filename,"%s/lx%d_ly%d_lz%d.tmp\0",dirname,lx,ly,lz);
	//printf("READING: L=(%d,%d,%d), filename=%s\n",lx,ly,lz,filename);
	sprintf(shellcommand,
		"if [ ! -e %s ]; then echo Reading Error: %s, does not exist; fi\0",
		filename,filename);
	fptr=fopen(filename,"r");
	fscanf(fptr,"%d %lf",&NRADIALread,&RADSTEPread);
	if(fabs(RADSTEPread-RADSTEP)>1.0E-6 || NRADIALread!=NRADIAL){
	  printf("Mesh in %s out of whack: RADSTEP=%g =? %g, NRADIAL=%d=?%d\n",
		 filename,RADSTEP,RADSTEPread,NRADIAL,NRADIALread);
	  exit(1);
	}
	for(ir=0;ir<NRADIAL;ir++){
	  fscanf(fptr,"%lf\n",&aa);
	  A[lx][ly][lz][ir]=aa;
	}
	fclose(fptr);
      }
    }
  }
  FillRemainderX();
  //printf("Reading finished\n");
}

void CCHArray::ReadAllA(char *dirname){
  char filename[160],dummy[200],shellcommand[200];
  int ir,lx,ly,lz,NRADIALread;
  double aa,RADSTEPread;
  FILE *fptr;

  for(lx=0;lx<=LMAX;lx+=dlx){
    for(ly=0;ly<=LMAX-lx;ly+=dly){
      for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
	sprintf(filename,"%s/lx%d_ly%d_lz%d.tmp\0",dirname,lx,ly,lz);
	//printf("READING: L=(%d,%d,%d), filename=%s\n",lx,ly,lz,filename);
	sprintf(shellcommand,
		"if [ ! -e %s ]; then echo Reading Error: %s, does not exist; fi\0",
		filename,filename);
	fptr=fopen(filename,"r");
	fscanf(fptr,"%d %lf",&NRADIALread,&RADSTEPread);
	if(fabs(RADSTEPread-RADSTEP)>1.0E-6 || NRADIALread!=NRADIAL){
	  printf("Mesh in %s out of whack: RADSTEP=%g =? %g, NRADIAL=%d=?%d\n",
		 filename,RADSTEP,RADSTEPread,NRADIAL,NRADIALread);
	  exit(1);
	}
	for(ir=0;ir<NRADIAL;ir++){
	  fscanf(fptr,"%lf\n",&aa);
	  A[lx][ly][lz][ir]=aa;
	}
	fclose(fptr);
      }
    }
  }
  FillRemainderX();
  //printf("Reading finished\n");
}

void CCHArray::PrintPars(){
  printf("XSYM=%d, YSYM=%d, ZSYM=%d\n",XSYM,YSYM,ZSYM);
  printf("LMAX=%d, NRADIAL=%d, RADSTEP=%g\n",LMAX,NRADIAL,RADSTEP);
}

void CCHArray::IncrementAExpArray(double x,double y,double z,
				       double weight){
  double ex,ey,ez;
  double r=sqrt(x*x+y*y+z*z);
  int ir=int(floor(r/RADSTEP));
  if(ir<NRADIAL){
    ex=x/r; ey=y/r; ez=z/r;
    IncrementAExpArrayFromE(ex,ey,ez,weight,ir);
  }
}

void CCHArray::IncrementAExpArrayFromE(double ex,double ey,double ez,
					double weight,int ir){
  int L,lx,ly,lz;
  double lfact;
  for(lx=0;lx<=1;lx+=dlx){
    for(ly=0;ly<=LMAX-lx;ly+=dly){
      for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
	L=lx+ly+lz;
	lfact=chcalc->DoubleFactorial(2*L+1)/chcalc->Factorial(L);
	IncrementElement(lx,ly,lz,ir,
			 weight*lfact*chcalc->GetAFromE(lx,ly,lz,ex,ey,ez));
      }
    }
  }
  FillRemainderX(ir);
}

void CCHArray::AltIncrementAExpArrayFromE(double ex,double ey,double ez,
						 double weight,int ir){
  int L,lx,ly,lz;
  double cL,cLL,delB,lfact;
  double ***B;
  B=new double**[2];
  for(lx=0;lx<2;lx++){
    B[lx]=new double *[LMAX+1-lx];
    for(ly=0;ly<=LMAX-lx;ly++){
      B[lx][ly]=new double[LMAX+1-lx-ly];
    }
  }
    
  B[0][0][0]=1.0;
  IncrementElement(0,0,0,ir,weight);
  for(L=1;L<=LMAX;L++){
    lfact=chcalc->DoubleFactorial(2*L+1)/chcalc->Factorial(L);
    cL=1.0/double(L);
    cLL=1.0/(L*(2*L-1));
    for(lx=0;lx<2;lx++){
      for(ly=0;ly<=L-lx;ly++){
	lz=L-lx-ly;
	delB=0.0;
	if(lx>0){
	  delB+=ex*(cL*lx-cLL*lx*(lx-1))*B[lx-1][ly][lz];
	}
	if(ly>0)
	  delB+=ey*(cL*ly-cLL*ly*(ly-1))*B[lx][ly-1][lz];
	if(lz>0)
	  delB+=ez*(cL*lz-cLL*lz*(lz-1))*B[lx][ly][lz-1];
	if(lx>1){
	  delB-=ey*cLL*lx*(lx-1)*B[lx-2][ly+1][lz];
	  delB-=ez*cLL*lx*(lx-1)*B[lx-2][ly][lz+1,0];
	}
	if(ly>1){
	  if(lx>0){
	    delB+=ex*cLL*ly*(ly-1)*(B[lx-1][ly-2][lz+2]+B[lx-1][ly][lz]);
	  }
	  else delB-=ex*cLL*ly*(ly-1)*B[lx+1][ly-2][lz];
	  delB-=ez*cLL*ly*(ly-1)*B[lx][ly-2][lz+1];
	}
	if(lz>1){
	  if(lx>0){
	    delB+=ex*cLL*lz*(lz-1)*(B[lx-1][ly+2][lz-2]+B[lx-1][ly][lz,0]);
	  }
	  else delB-=ex*cLL*lz*(lz-1)*B[lx+1][ly][lz-2];
	  delB-=ey*cLL*lz*(lz-1)*B[lx][ly+1][lz-2];
	}
	B[lx][ly][lz]=delB;
	if(lx%dlx==0 && ly%dly==0 && lz%dlz==0){
	  IncrementElement(lx,ly,lz,ir,delB*lfact*weight);
	}
      }
    }
    //FillRemainderX(B,0);
  }
  for(lx=0;lx<2;lx++){
    for(ly=0;ly<=LMAX-lx;ly++){
      delete [] B[lx][ly];
    }
    delete [] B[lx];
  }
  delete [] B;
}

double CCHArray::GetAExpElementFromMArray(int lx,int ly,int lz,int ir){
  int L,m,mx,my,mz;
  double factor,answer,lfact;
  L=lx+ly+lz;
  lfact=chcalc->DoubleFactorial(2*L+1)/chcalc->Factorial(L);
  answer=0.0;
  for(mx=0;mx<=lx/2;mx++){
    for(my=0;my<=ly/2;my++){
      for(mz=0;mz<=lz/2;mz++){
	m=mx+my+mz;
	if(m>0)	factor=pow(-0.5,m);
	else factor=1.0;
	if(L>m) factor*=chcalc->DoubleFactorial(2*L-2*m-1);
	if(L>0) factor=factor/chcalc->DoubleFactorial(2*L-1);
	factor*=chcalc->Factorial(lx)/(chcalc->Factorial(lx-2*mx)*chcalc->Factorial(mx));
	factor*=chcalc->Factorial(ly)/(chcalc->Factorial(ly-2*my)*chcalc->Factorial(my));
	factor*=chcalc->Factorial(lz)/(chcalc->Factorial(lz-2*mz)*chcalc->Factorial(mz));
	answer+=factor*lfact*GetElement(lx-2*mx,ly-2*my,lz-2*mz,ir);
      }
    }
  }
  return answer;
}

double CCHArray::GetMElementFromAExpArray(int lx,int ly,int lz,int ir){
  int mx,my,mz,m,L;
  double factor,answer,*lfact;
  L=lx+ly+lz;
  lfact=new double[L+1];
  lfact[0]=1.0;
  for(m=1;m<=L;m++) lfact[m]=lfact[m-1]*double(m)/double(2*m+1);
  answer=0.0;
  for(mx=0;mx<=lx/2;mx++){
    for(my=0;my<=ly/2;my++){
      for(mz=0;mz<=lz/2;mz++){
	m=mx+my+mz;
	factor=pow(0.5,m);
	factor=factor*chcalc->DoubleFactorial(2*L-4*m+1)/chcalc->DoubleFactorial(2*L-2*m+1);
	factor*=chcalc->Factorial(lx)/(chcalc->Factorial(lx-2*mx)*chcalc->Factorial(mx));
	factor*=chcalc->Factorial(ly)/(chcalc->Factorial(ly-2*my)*chcalc->Factorial(my));
	factor*=chcalc->Factorial(lz)/(chcalc->Factorial(lz-2*mz)*chcalc->Factorial(mz));
	answer+=factor*lfact[L-2*m]
	  *GetElement(lx-2*mx,ly-2*my,lz-2*mz,ir);
      }
    }
  }
  delete lfact;
  return answer;
}	     

void CCHArray::IncrementAExpArrayFromThetaPhi(double theta,double phi,
					       double weight,int ir){
  double stheta,ex,ey,ez;
  stheta=sin(theta);
  ex=stheta*cos(phi);
  ey=stheta*sin(phi);
  ez=cos(theta);
  IncrementAExpArrayFromE(ex,ey,ez,weight,ir);
}

void CCHArray::IncrementMArrayFromE(double ex,double ey,double ez,double weight,
				    int ir){
  int lx,ly,lz;
  for(lx=0;lx<=LMAX;lx+=dlx){
    for(ly=0;ly<=LMAX-lx;ly+=dly){
      for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
	IncrementElement(lx,ly,lz,ir,pow(ex,lx)*pow(ey,ly)*pow(ez,lz)*weight);
      }
    }
  }
}

void CCHArray::IncrementMArrayFromThetaPhi(double theta,double phi,
					   double weight,int ir){
  double stheta,ex,ey,ez;
  stheta=sin(theta);
  ex=stheta*cos(phi);
  ey=stheta*sin(phi);
  ez=cos(theta);
  IncrementMArrayFromE(ex,ey,ez,weight,ir);
}

void CCHArray::FillRemainderX(){
  for(int ir=0;ir<NRADIAL;ir++) FillRemainderX(ir);
}
void CCHArray::FillRemainderY(){
  for(int ir=0;ir<NRADIAL;ir++) FillRemainderY(ir);
}
void CCHArray::FillRemainderZ(){
  for(int ir=0;ir<NRADIAL;ir++) FillRemainderZ(ir);
}

void CCHArray::FillRemainderX(int ir){
  int L,lx,ly,lz;
  if(LMAX>1){
    for(lx=2;lx<=LMAX;lx+=dlx){
      for(ly=0;ly<=LMAX-lx;ly+=dly){
	for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
	  SetElement(lx,ly,lz,ir,-GetElement(lx-2,ly+2,lz,ir)
			-GetElement(lx-2,ly,lz+2,ir));
	}
      }
    }
  }
}

void CCHArray::FillRemainderY(int ir){
  int lx,ly,lz;
  if(LMAX>1){
    for(ly=2;ly<=LMAX;ly+=dly){
      for(lx=0;lx<=LMAX-ly;lx+=dlx){
	for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
	  SetElement(lx,ly,lz,ir,-GetElement(lx,ly-2,lz+2,ir)
			-GetElement(lx+2,ly-2,lz,ir));
	}
      }
    }
  }
}

void CCHArray::FillRemainderZ(int ir){
  int lx,ly,lz;
  if(LMAX>1){
    for(lz=2;lz<=LMAX;lz+=dlz){
      for(lx=0;lx<=LMAX-lz;lx+=dlx){
	for(ly=0;ly<=LMAX-lx-lz;ly+=dly){
	  SetElement(lx,ly,lz,ir,-GetElement(lx,ly+2,lz-2,ir)
			-GetElement(lx+2,ly,lz-2,ir));

	}
      }
    }
  }
}

double CCHArray::AExpand(double theta,double phi,int ir){
  double ex,ey,ez,sthet;
  sthet=sin(theta);
  ex=sthet*cos(phi);
  ey=sthet*sin(phi);
  ez=cos(theta);
  return AExpand(ex,ey,ez,ir);
}

double CCHArray::AExpand(double ex,double ey,double ez,int ir){
  double dela,answer=0.0;
  int lx,ly,lz;

  for(lx=0;lx<=LMAX;lx+=dlx){
    for(ly=0;ly<=LMAX-lx;ly+=dly){
      for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
	dela=GetElement(lx,ly,lz,ir)*pow(ex,lx)*pow(ey,ly)*pow(ez,lz);
	answer+=dela*chcalc->Trinomial(lx,ly,lz);
      }
    }
  }
  return answer;
}
double CCHArray::AExpand(double x,double y,double z){
  double r,ex,ey,ez;
  int ir;
  double answer=0.0;
  r=sqrt(x*x+y*y+z*z);
  ir=int(floor(r/RADSTEP));
  if(ir<NRADIAL){
    ex=x/r; ey=y/r; ez=z/r;
    answer=AExpand(ex,ey,ez,ir);
  }
  return answer;
}


void CCHArray::RandomInit(int iseed){
  randy=new CRandom(iseed);
}

void CCHArray::Randomize(double mag,int ir){
  int lx,ly,lz;
  double value;
  if(randy==NULL) RandomInit(-12345);
  for(lx=0;lx<=LMAX;lx+=dlx){
    for(ly=0;ly<=LMAX-lx;ly+=dly){
      for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
	value=(1.0-2.0*randy->ran())*mag;
	SetElement(lx,ly,lz,ir,value);
      }
    }
  }  
}
void CCHArray::Randomize(double mag){
  int lx,ly,lz,ir;
  double value;
  if(randy==NULL) RandomInit(-12345);
  for(ir=0;ir<NRADIAL;ir++){
    Randomize(mag,ir);
  }  
}

void CCHArray::RandomizeA(double mag,int ir){
  int lx,ly,lz;
  double value;
  if(randy==NULL) RandomInit(-12345);
  for(lx=0;lx<=1;lx+=dlx){
    for(ly=0;ly<=LMAX-lx;ly+=dly){
      for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
	value=(1.0-2.0*randy->ran())*mag;
	SetElement(lx,ly,lz,ir,value);
      }
    }
    FillRemainderX(ir);
  }
}
void CCHArray::RandomizeA(double mag){
  int lx,ly,lz,ir;
  double value;
  if(randy==NULL) RandomInit(-12345);
  for(ir=0;ir<NRADIAL;ir++){
    RandomizeA(mag,ir);
  }
}

void CCHArray::RandomizeA_Gaussian(double mag,int ir){
  int lx,ly,lz;
  double value;
  if(randy==NULL) RandomInit(-12345);
  for(lx=0;lx<=1;lx+=dlx){
    for(ly=0;ly<=LMAX-lx;ly+=dly){
      for(lz=0;lz<=LMAX-lx-ly;lz+=dlz){
	value=randy->gauss()*mag;
	SetElement(lx,ly,lz,ir,value);
      }
    }
    FillRemainderX(ir);
  }
}
void CCHArray::RandomizeA_Gaussian(double mag){
  int lx,ly,lz,ir;
  double value;
  if(randy==NULL) RandomInit(-12345);
  for(ir=0;ir<NRADIAL;ir++){
    RandomizeA_Gaussian(mag,ir);
  }
}

void CCHArray::Detrace(int ir){
  CCHArray *B;
  B=new CCHArray(LMAX,1,RADSTEP,XSYM,YSYM,ZSYM);
  ArrayCalc::Detrace(this,ir,B,0);
  ArrayCalc::CopyArray(B,0,this,ir);
  delete(B);
}

void CCHArray::Detrace(){
  int ir;
  for(ir=0;ir<NRADIAL;ir++) Detrace(ir);
}

void CCHArray::WriteShort(char *filename,int LGMAX){
    int ir,lx,ly,lz,L;
  FILE *fptr;
  if(LGMAX>LMAX) LGMAX=LMAX;
  fptr=fopen(filename,"w");
  fprintf(fptr,"!       ");
  for(L=0;L<=LGMAX;L++){
    for(lx=0;lx<=1;lx++){
      for(ly=0;ly<=L-lx;ly++){
	lz=L-lx-ly;
	if(lx%dlx==0 && ly%dly==0 && lz%dlz==0)
	  fprintf(fptr,"  (%d,%d,%d)  ",lx,ly,lz);
      }
    }
  }
  fprintf(fptr,"\n");
  for(ir=0;ir<NRADIAL;ir++){
    fprintf(fptr,"%7.2f ",RADSTEP*(ir+0.5));
    for(L=0;L<=LGMAX;L++){
      for(lx=0;lx<=1;lx++){
	for(ly=0;ly<=L-lx;ly++){
	  lz=L-lx-ly;
	  if(lx%dlx==0 && ly%dly==0 && lz%dlz==0)
	    fprintf(fptr," %9.2e ",GetElement(lx,ly,lz,ir));
	}
      }
    }
    fprintf(fptr,"\n");
  }
  fclose(fptr);

}

#endif