#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.14159265358979323844
#define NB 1

void gauss(double *g1,double *g2);
long IDUM=-1234;double ran2(void);
#define NTYPES 5

int main()
{
  double p[4],r[4],mass[NTYPES]={938.3,939.6,139.58,139.58,493.7};
  int ident[NTYPES]={2212,2112,211,-211,321};
  double dummy,weight,pmag2,pmag,temp,rx,ry,rz,tau,cthet,sthet,phi;
  int i,id,alpha,iblankline,ib;
  char filename[80];
  FILE *fptr;
  printf("idents are %d %d %d\n",ident[0],ident[1],ident[2]);
  printf("What is the temperature? (in MeV)\n");
  scanf("%lf",&temp);
  printf("What are rx,ry,rz and tau? (in fm)\n");
  scanf("%lf %lf %lf %lf",&rx,&ry,&rz,&tau);
  printf("radii are %g %g %g %g\n",rx,ry,rz,tau);
  for(ib=0;ib<NB;ib++)
    {
      sprintf(filename,"thermal_gauss%02d.dat",ib+1);
      fptr=fopen(filename,"w");
      for(iblankline=0;iblankline<3;iblankline++)
        {
          fprintf(fptr,"line number %d, blah blah blah\n",iblankline);
        }
      fprintf(fptr,"%d %d %lf %lf\n",1,10000,0.0,0.0);
      for(i=0;i<10000;i++)
        {
          id=(int)floor(NTYPES*ran2());
          if(mass[id]/temp>8.0)
            {
              pmag2=0.0;
              gauss(&p[1],&p[2]);gauss(&p[3],&dummy);
              for(alpha=1;alpha<4;alpha++)
                {
                  p[alpha]=sqrt(temp*mass[id])*p[alpha];
                  pmag2=pmag2+p[alpha]*p[alpha];
                }
              p[0]=sqrt(pmag2+mass[id]*mass[id]);
            }
          else
            {
              do
                {
                  pmag=-temp*log(ran2()*ran2()*ran2());
                  p[0]=sqrt(mass[id]*mass[id]+pmag*pmag);
                  weight=exp(-(p[0]-pmag)/temp);
                }while(ran2()>weight);
              //printf("id=%d  pmag=%g\n",ident[id],pmag);
              cthet=1.0-2.0*ran2();
              sthet=sqrt(1.0-cthet*cthet);
              phi=2.0*pi*ran2();
              p[3]=pmag*cthet;
              p[1]=pmag*sthet*cos(phi);
              p[2]=pmag*sthet*sin(phi);
            }
          gauss(&r[0],&r[1]);gauss(&r[2],&r[3]);
          r[0]=tau*r[0];
          r[1]=rx*r[1];
          r[2]=ry*r[2];
          r[3]=rz*r[3];
          fprintf(fptr,"%d %d %9g %9g %9g %9g %9g %9g %9g %9g %9g\n",
                  i+1,ident[id],p[1]/1000.0,p[2]/1000.0,p[3]/1000.0,p[0]/1000.0,
                  mass[id]/1000.0,r[1],r[2],r[3],r[0]);
        }
      fclose(fptr);
    }
  return 0;
}

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2()
{
  long J;
  long K;
  static long IDUM2=123456789;
  static long IY=0;
  static long IV[NTAB];
  double randy;
  if(IDUM <=0)
    {
      if(-(IDUM)<1) IDUM=1;
      else IDUM = -(IDUM);
      IDUM2=(IDUM);
      for (J=NTAB+7;J>=0;J--)
        {
          K=(IDUM)/IQ1;
          IDUM=IA1*(IDUM-K*IQ1)-K*IR1;
          if(IDUM<0) IDUM +=IM1;
          if(J<NTAB) IV[J]=IDUM;
        }
      IY=IV[0];
      //printf("Initialization Completed %d\n",IDUM);
    }
  K=(IDUM)/IQ1;
  IDUM=IA1*(IDUM-K*IQ1)-K*IR1;
  if(IDUM<0) IDUM+=IM1;
  K=IDUM2/IQ2;
  IDUM2=IA2*(IDUM2-K*IQ2)-K*IR2;
  if(IDUM2<0) IDUM2+=IM2;
  J=IY/NDIV;
  IY=IV[J]-IDUM2;
  IV[J]=IDUM;
  if(IY<1) IY+=IMM1;
  randy=AM*IY;
  if(randy>RNMX) return RNMX;
  else return randy;
}
//**********************************************************
void gauss(double *g1,double *g2)
{
  double x,y,r2,r;
  do
    {
      x=1.0-2.0*ran2();
      y=1.0-2.0*ran2();
      r2=x*x+y*y;
    }
  while (r2>1.0);
  r=sqrt(r2);
  *g1=(x/r)*sqrt(-2.0*log(r2));
  *g2=(y/x)**g1;
}

