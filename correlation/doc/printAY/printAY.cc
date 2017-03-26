#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include "arrays.cc"
#include "parametermap.cc"
#include "sf.cc"
#include "misc.cc"
#include "random.cc"

using namespace std;
void PrintAofY_tex();
void PrintYofA_tex();

bool remainder(double x, int i);
void cfactorize(int L,double &crealnum,double &cimagnum,double &crealdenom,
		double &cimagdenom,double &csqrtnum,double  &csqrtdenom);
int doublefact(int);
int fact(int);

int main(){
  PrintAofY_tex();
  PrintYofA_tex();
  return 0;
}
void PrintYofA_tex(){
  int PrintLmax=4;
  int L,lx,ly,lz,m,k,kk,modd;
  complex<double> c;
  double csqrtnum,csqrtdenom;
  double crealnum,cimagnum,crealdenom,cimagdenom;
  bool firstterm,success;
  double csqrt;
  FILE *fptr;
  complex<double> ci(0.0,1.0); 
  fptr=fopen("YofA.tex","w");
  fprintf(fptr,"\\documentclass{article}\n");
  fprintf(fptr,"\\parindent 0pt\n");
  fprintf(fptr,"\\begin{document}\n");
  for(L=1;L<=PrintLmax;L++){
    for(m=0;m<=L;m++){
      csqrt=double(2*L+1)/double(fact(L+m)*fact(L-m));
      if(m>0)fprintf(fptr,"$Y_{%d,\\pm %d} =",L,m);
      else fprintf(fptr,"$Y_{%d,0} =",L);
      firstterm=0;
      for(k=0;k<=m;k++){ 
	c=pow(-1.0,m)*pow(ci,m-k)*double(doublefact(2*L-1))
	  *(double(fact(m))/double(fact(k)*fact(m-k)));
	csqrtnum=(2*L+1);
	csqrtdenom=4*fact(L+m)*fact(L-m);
	crealnum=real(c);
	cimagnum=imag(c);
	crealdenom=1;
	cimagdenom=1;

	cfactorize(L,crealnum,cimagnum,crealdenom,cimagdenom,
		   csqrtnum,csqrtdenom);
	modd=m%2;
	if(fabs(crealnum*crealnum+cimagnum*cimagnum)>1.0E-8){
	  fprintf(fptr,"\n");
	  
	  if(fabs(crealnum)>1.0E-8){
	    if(modd==0 && firstterm !=0 && crealnum>1.0E-8) 
	      fprintf(fptr," +");
	    if(modd==0 && crealnum<-1.0E-8) fprintf(fptr," -");
	    if(modd==1 && crealnum>1.0E-8) fprintf(fptr," \\pm");
	    if(modd==1 && crealnum<-1.0E-8) fprintf(fptr," \\mp");
	    if(fabs(fabs(crealdenom)-1.0)>1.0E-8) fprintf(fptr,"(");
	    if(fabs(crealnum)-1>1.0E-8) fprintf(fptr,"%g",fabs(crealnum));
	    else if(fabs(fabs(crealdenom)-1)>1.0E-8) fprintf(fptr,"1");
	    if(fabs(fabs(crealdenom)-1.0)>1.0E-8) 
	      fprintf(fptr,"/%g)",crealdenom);
	  }
	  if(fabs(cimagnum)>1.0E-8){
	    if(modd==1 && firstterm !=0 && cimagnum>1.0E-8)
	      fprintf(fptr," +");
	    if(modd==1 && cimagnum<-1.0E-8) fprintf(fptr," -");
	    if(modd==0 && cimagnum>1.0E-8) fprintf(fptr," \\pm");
	    if(modd==0 && cimagnum<-1.0E-8) fprintf(fptr," \\mp");
	    if(fabs(fabs(cimagdenom)-1.0)>1.0E-8) fprintf(fptr,"(");
	    if(fabs(fabs(cimagnum)-1.0)<1.0E-8) fprintf(fptr,"i");
	    else fprintf(fptr,"%gi",fabs(cimagnum));
	    if(fabs(fabs(cimagdenom)-1.0)>1.0E-8)
	      fprintf(fptr,"/%g)",cimagdenom);
	  }
	  fprintf(fptr,"\\sqrt{");
	  fprintf(fptr,"%g",csqrtnum);
	  if(fabs(fabs(csqrtdenom)-1.0)>1.0E-8) 
	    fprintf(fptr,"/%g\\pi}",csqrtdenom);
	  else fprintf(fptr,"/\\pi}");
	  fprintf(fptr,"{\\mathcal A}^{(%d)}_{",L);
	  for(kk=0;kk<k;kk++) fprintf(fptr,"x");
	  for(kk=0;kk<m-k;kk++) fprintf(fptr,"y");
	  for(kk=0;kk<L-m;kk++) fprintf(fptr,"z");
	  fprintf(fptr,"}");
	  firstterm=1;
	}
      }
      fprintf(fptr,"$\\\\\n");
    }

  }
  fprintf(fptr,"\\end{document}\n");
  fclose(fptr);
}

bool remainder(double x,int i){
  if(fabs(x)>1.0E-8 && fabs(x-i*rint(x/double(i)))<1.0E-8) return 1;
  else return 0;
}

void cfactorize(int L,double &crealnum,double &cimagnum,double &crealdenom,
		double &cimagdenom,double &csqrtnum,double  &csqrtdenom){
  int k,success;
  double EPSILON=1.0E-8;
  do{
    success=0;
    for(k=2;k<=2*L+1;k++){
      if(remainder(csqrtnum,k) && remainder(csqrtdenom,k)){
	csqrtnum=csqrtnum/double(k);
	csqrtdenom=csqrtdenom/double(k);
	success=1;
      }
    }
  
    for(k=2;k<=2*L+1;k++){
      if(fabs(csqrtnum)>EPSILON && remainder(csqrtnum,k*k)){
	csqrtnum=csqrtnum/double(k*k);
	crealnum=crealnum*k;
	cimagnum=cimagnum*k;
	success=1;
      }
      if(fabs(csqrtdenom)>EPSILON && remainder(csqrtdenom,k*k)){
	csqrtdenom=csqrtdenom/double(k*k);
	crealdenom=crealdenom*k;
	cimagdenom=cimagdenom*k;
	success=1;
      }
    }
  
    for(k=2;k<=2*L+1;k++){
      if(remainder(crealnum,k) && remainder(crealdenom,k)){
	crealnum=crealnum/double(k);
	crealdenom=crealdenom/double(k);
	success=1;
      }
      if(remainder(cimagnum,k) && remainder(cimagdenom,k)){
	cimagnum=cimagnum/double(k);
	cimagdenom=cimagdenom/double(k);
	success=1;
      }
    }
    for(k=2;k<=2*L+1;k++){
      if(fabs(crealnum)>EPSILON){
	if(remainder(csqrtnum,k) && remainder(crealdenom,k)){
	  csqrtnum=csqrtnum/double(k);
	  csqrtdenom=csqrtdenom*double(k);
	  crealdenom=crealdenom/double(k);
	  success=1;
	}
      }
      if(fabs(cimagnum)>EPSILON){
	if(remainder(csqrtnum,k) && remainder(cimagdenom,k)){
	  csqrtnum=csqrtnum/double(k);
	  csqrtdenom=csqrtdenom*double(k);
	  cimagdenom=cimagdenom/double(k);
	  success=1;
	}
      }
    }
    for(k=2;k<=2*L+1;k++){
      if(fabs(crealnum)>EPSILON){
	if(remainder(csqrtdenom,k) && remainder(crealnum,k)){
	  csqrtdenom=csqrtdenom/double(k);
	  csqrtnum=csqrtnum*double(k);
	  crealnum=crealnum/double(k);
	  success=1;
	}
      }
      if(fabs(cimagnum)>EPSILON){
	if(remainder(csqrtdenom,k) && remainder(cimagnum,k)){
	  csqrtdenom=csqrtdenom/double(k);
	  csqrtnum=csqrtnum*double(k);
	  cimagnum=cimagnum/double(k);
	  success=1;
	}
      }
    }
  }while(success==1);	      
}

void PrintAofY_tex(){
  CCHCalc *chcalc=new CCHCalc();
  const double PI=4.0*atan(1.0);
  int PrintLmax=4;
  int L,lx,ly,lz,m,k,kk;
  complex<double> *c;
  double csqrtnum,csqrtdenom;
  double crealnum,cimagnum,crealdenom,cimagdenom;
  bool firstterm,success;
  double *csqrt;
  FILE *fptr;
  complex<double> ci(0.0,1.0); 
  fptr=fopen("AofY.tex","w");
  fprintf(fptr,"\\documentclass{article}\n");
  fprintf(fptr,"\\parindent 0pt\n");
  fprintf(fptr,"\\begin{document}\n");
  c=new complex<double>[PrintLmax+1];
  csqrt=new double[PrintLmax+1];
  for(L=1;L<=PrintLmax;L++){
    for(lx=0;lx<=L;lx++){
      for(ly=0;ly<=L-lx;ly++){
	lz=L-lx-ly;
	for(m=L;m>=0;m--){
	  csqrt[m]=double(2*L+1)/double(fact(L+m)*fact(L-m));
	  c[m]=0.0;
	  for(k=0;k<=m;k++) 
	    c[m]+=pow(-ci,m-k)*(double(fact(m))/double(fact(k)*fact(m-k)))
	      *chcalc->GetOverlap0(lx,ly,lz,k,m-k,L-m)/(4.0*PI);
	  c[m]=c[m]*pow(-1.0,m)*double(doublefact(2*L-1));
	}
	fprintf(fptr,"${\\mathcal A}^{(%d)}_{",L);
	for(kk=0;kk<lx;kk++) fprintf(fptr,"x");
	for(kk=0;kk<ly;kk++) fprintf(fptr,"y");
	for(kk=0;kk<lz;kk++) fprintf(fptr,"z");
	fprintf(fptr,"} =");
	firstterm=0;
	for(m=L;m>=-L;m--){
	  csqrtnum=4*csqrt[abs(m)]*fact(2*L+1);
	  csqrtdenom=fact(2*L+1);
	  crealnum=real(c[abs(m)])*fact(2*L+1);
	  cimagnum=imag(c[abs(m)])*fact(2*L+1);
	  crealdenom=fact(2*L+1);
	  cimagdenom=fact(2*L+1);
	  if(m<0){
	    crealnum=crealnum*pow(-1.0,-m);
	    cimagnum=cimagnum*pow(-1.0,-m);
	    if(abs(ly)%2!=0){
	      crealnum*=-1.0;
	      cimagnum*=-1.0;
	    }
	  }

	  cfactorize(L,crealnum,cimagnum,crealdenom,cimagdenom,
		     csqrtnum,csqrtdenom);

	  if(fabs(crealnum*crealnum+cimagnum*cimagnum)>1.0E-8){
	    fprintf(fptr,"\n");
	    if(m<L && firstterm==1 && 
	       (crealnum>1.0E-8 || cimagnum>1.0E-8)){
	      fprintf(fptr," +");
	    }
	    if(crealnum<-1.0E-8 || cimagnum<-1.0E-8) fprintf(fptr," -");
	    firstterm=1;

	    if(fabs(crealnum)>1.0E-8){
	      if(fabs(fabs(crealdenom)-1.0)>1.0E-8) fprintf(fptr,"(");
	      if(fabs(crealnum)-1>1.0E-8) fprintf(fptr,"%g",fabs(crealnum));
	      else if(fabs(fabs(crealdenom)-1)>1.0E-8) fprintf(fptr,"1");
	      if(fabs(fabs(crealdenom)-1.0)>1.0E-8) 
		fprintf(fptr,"/%g)",crealdenom);
	    }
	    if(fabs(cimagnum)>1.0E-8){
	      if(fabs(fabs(cimagdenom)-1.0)>1.0E-8) fprintf(fptr,"(");
	      if(fabs(fabs(cimagnum)-1.0)<1.0E-8) fprintf(fptr,"i");
	      else fprintf(fptr,"%gi",fabs(cimagnum));
	      if(fabs(fabs(cimagdenom)-1.0)>1.0E-8)
		fprintf(fptr,"/%g)",cimagdenom);
	    }
	    fprintf(fptr,"\\sqrt{");
	    if(fabs(csqrtnum-1)>1.0E-8) fprintf(fptr,"%g",csqrtnum);
	    fprintf(fptr,"\\pi");
	    if(csqrtdenom>1.0001) fprintf(fptr,"/%g",csqrtdenom);
	    fprintf(fptr,"}Y_{%d,%d}",L,m);
	  }
	}
	fprintf(fptr,"$\\\\\n");
      }
    }
    //fprintf(fptr,"______________________________________________\n");
  }
  fprintf(fptr,"\\end{document}\n");
  fclose(fptr);
}


int doublefact(int n){
  int fact=1;
  int i;
  if(n>1)
    for(i=n;i>=1;i-=2) fact*=i;
  else
    fact=1;
  return fact;
}
int fact(int n){
  int fact=1;
  int i;
  if(n>1)
    for(i=n;i>=1;i-=1) fact*=i;
  else
    fact=1;
  return fact;
}


