#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
using namespace std;
#include "sf.cc"
#include "misc.cc"


int main(){
  const double PI=4.0*atan(1.0);
  int nqmax=700;
  int iq,iread,iqsmooth,ichannel;
  double qsmooth,q,a,elab,plab,roots,eta,deltaread[61],qread[61];
  double delta[3][nqmax+1],ddeltadq[3][nqmax+1];
  char dummy[200];
  FILE *fptr;
  int ell[3]={0,1,1};
  double e1,e2,e,MPROTON=938.28,MKAON=493.677;
  double delq=0.5;

  iqsmooth=1;
  fptr=fopen("s11.dat\0","r");
  fgets(dummy,200,fptr);
  for(iread=0;iread<=60;iread++){
    fscanf(fptr,"%lf %lf",&plab,&deltaread[iread]);
    fgets(dummy,200,fptr);
    if(iread>0){
      elab=MPROTON+sqrt(MKAON*MKAON+plab*plab);
      roots=sqrt(elab*elab-plab*plab);
      qread[iread]=sqrt(Misc::triangle(roots,MPROTON,MKAON));
    }
    else qread[iread]=0.0;
    deltaread[iread]=deltaread[iread]*PI/180.0;
    if(iread==iqsmooth){
      a=tan(deltaread[iread])/qread[iread];
      qsmooth=qread[iread];
    }
  }
  fclose(fptr);
  
  iread=0;
  for(iq=0;iq<=nqmax;iq++){
    q=delq*double(iq);;
    if(q>qsmooth){
      if(iread>0) iread=iread-1;
      do{
	iread+=1;
      }while(q>qread[iread] && iread<=60);
      delta[0][iq]=((q-qread[iread-1])*deltaread[iread]
		    +((qread[iread]-q)*deltaread[iread-1]))
	/(qread[iread]-qread[iread-1]);
      ddeltadq[0][iq]=(deltaread[iread]-deltaread[iread-1])
	/(qread[iread]-qread[iread-1]);
    }
    else{
      delta[0][iq]=atan(fabs(a*q))*a/fabs(a);
      ddeltadq[0][iq]=a*pow(cos(delta[0][iq]),2);
    }
  }

  // ________________________________________________________

  iqsmooth=5;
  fptr=fopen("p11.dat\0","r");
  fgets(dummy,200,fptr);
  for(iread=0;iread<=60;iread++){
    fscanf(fptr,"%lf %lf",&plab,&deltaread[iread]);
    fgets(dummy,200,fptr);
    if(iread>0){
      elab=MPROTON+sqrt(MKAON*MKAON+plab*plab);
      roots=sqrt(elab*elab-plab*plab);
      qread[iread]=sqrt(Misc::triangle(roots,MPROTON,MKAON));
    }
    else qread[iread]=0.0;
    deltaread[iread]=deltaread[iread]*PI/180.0;
    if(iread==iqsmooth){
      a=tan(deltaread[iread])/pow(qread[iread],3);
      qsmooth=qread[iread];;
    }
  }
  fclose(fptr);

  iread=0;
  for(iq=0;iq<=nqmax;iq++){
    q=delq*double(iq);
    if(q>qsmooth){
      if(iread>0) iread=iread-1;
      do{
	iread+=1;
      }while(q>qread[iread] && iread<=60);
      delta[1][iq]=((q-qread[iread-1])*deltaread[iread]
		    +((qread[iread]-q)*deltaread[iread-1]))
	/(qread[iread]-qread[iread-1]);
      ddeltadq[1][iq]=(deltaread[iread]-deltaread[iread-1])
	/(qread[iread]-qread[iread-1]);
    }
    else{
      delta[1][iq]=atan(fabs(a*q*q*q))*a/fabs(a);
      ddeltadq[1][iq]=3*a*q*q*pow(cos(delta[1][iq]),2);
    }
  }

  // ________________________________________________________

  iqsmooth=5;
  fptr=fopen("p13.dat\0","r");
  fgets(dummy,200,fptr);
  for(iread=0;iread<=60;iread++){
    fscanf(fptr,"%lf %lf",&plab,&deltaread[iread]);
    fgets(dummy,200,fptr);
    if(iread>0){
      elab=MPROTON+sqrt(MKAON*MKAON+plab*plab);
      roots=sqrt(elab*elab-plab*plab);
      qread[iread]=sqrt(Misc::triangle(roots,MPROTON,MKAON));
      //printf("iread=%d, qread=%g, delta=%g\n",iread,qread[iread],deltaread[iread]);
    }
    else qread[iread]=0.0;
    deltaread[iread]=deltaread[iread]*PI/180.0;
    if(iread==iqsmooth){
      a=tan(deltaread[iread])/pow(qread[iread],3);
      qsmooth=qread[iread];
    }
  }
  fclose(fptr);

  iread=0;
  for(iq=0;iq<=nqmax;iq++){
    q=double(iq)*delq;
    if(q>qsmooth){
      if(iread>0) iread=iread-1;
      do{
	iread+=1;
      }while(q>qread[iread] && iread<=60);
      delta[2][iq]=((q-qread[iread-1])*deltaread[iread]
		    +((qread[iread]-q)*deltaread[iread-1]))
	/(qread[iread]-qread[iread-1]);
      ddeltadq[2][iq]=(deltaread[iread]-deltaread[iread-1])
	/(qread[iread]-qread[iread-1]);
    }
    else{
      delta[2][iq]=atan(fabs(a*q*q*q))*a/fabs(a);
      ddeltadq[2][iq]=3*a*q*q*pow(cos(delta[2][iq]),2);
    }
  }

  for(iq=1;iq<=nqmax;iq++){
    q=delq*double(iq);
    e1=sqrt(q*q+MPROTON*MPROTON);
    e2=sqrt(q*q+MKAON*MKAON);
    e=e1+e2;
    //eta=e1*e2/((e1+e2)*137.036*q);// old way is wrong
    eta=(pow(e,4)-pow(MPROTON*MPROTON-MKAON*MKAON,2))/(4.0*e*e*e*137.036*q);

    for(ichannel=0;ichannel<3;ichannel++)
      CoulWave::phaseshift_CoulombCorrect(ell[ichannel],q,eta,
					  delta[ichannel][iq],
					  ddeltadq[ichannel][iq]);
  }
  for(ichannel=0;ichannel<3;ichannel++){
    delta[ichannel][0]=ddeltadq[ichannel][0]=0.0;
  }

  fptr=fopen("pk_phaseshiftdat.cc","w");
  fprintf(fptr,"   double delqdata=%g;\n",delq);
  fprintf(fptr,"   int nqdata=%d;\n",nqmax);
  fprintf(fptr,"   double data_delta[%d][%d]={",3,nqmax+1);
  for(ichannel=0;ichannel<3;ichannel++){
    fprintf(fptr,"{0");
    for(iq=1;iq<=nqmax;iq++) fprintf(fptr,",%g",delta[ichannel][iq]);
    fprintf(fptr,"}");
    if(ichannel!=2) fprintf(fptr,",\n      ");
  }
  fprintf(fptr,"};\n");
  fprintf(fptr,"   double data_ddeltadq[%d][%d]={",3,nqmax+1);
  for(ichannel=0;ichannel<3;ichannel++){
    fprintf(fptr,"{%g",ddeltadq[0][0]);
    for(iq=1;iq<=nqmax;iq++) fprintf(fptr,",%g",ddeltadq[ichannel][iq]);
    fprintf(fptr,"}");
    if(ichannel!=2) fprintf(fptr,",\n      ");
  }
  fprintf(fptr,"};\n");
  fclose(fptr);

  fptr=fopen("pk_phaseshifts.dat","w");
  for(iq=0;iq<nqmax;iq++){
    q=delq*iq;
    fprintf(fptr,"%5.2f ",q);
    for(ichannel=0;ichannel<3;ichannel++){
      fprintf(fptr,"%8.5lf %8.5lf ",delta[ichannel][iq],ddeltadq[ichannel][iq]);
    }
    fprintf(fptr,"\n");
  }
  fclose(fptr);

}
