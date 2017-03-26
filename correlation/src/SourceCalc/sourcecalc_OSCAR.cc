#ifndef __INCLUDE_SOURCECALC_OSCAR
#define __INCLUDE_SOURCECALC_OSCAR

using namespace std;

CSourceCalc_OSCAR::CSourceCalc_OSCAR(){
  InitSPars();
}

void CSourceCalc_OSCAR::InitSPars(){
  // DEFAULT VALUES
  parameter::set(spars,"Pt",600);
  parameter::set(spars,"DELPT",20);
  parameter::set(spars,"PHIMIN_DEG",0);
  parameter::set(spars,"PHIMAX_DEG",360.0);
  parameter::set(spars,"YMIN",-3.0);
  parameter::set(spars,"YMAX",3.0);
  parameter::set(spars,"IDa",211);
  parameter::set(spars,"IDb",211);
  parameter::set(spars,"Ma",139.58);
  parameter::set(spars,"Mb",139.58);
  parameter::set(spars,"AEQUALB",0);
  parameter::set(spars,"NMAX",20000);
  parameter::set(spars,"OSCARfilename","UNDEFINED");
  parameter::set(spars,"NEVENTSMAX",10000);
}

void CSourceCalc_OSCAR::SetSPars(double Pt_set,double DELPT_set,
			     double PHIMIN_DEG_set,double PHIMAX_DEG_set,
			     double YMIN_set,double YMAX_set){
  parameter::set(spars,"Pt",Pt_set);
  parameter::set(spars,"DELPT",DELPT_set);
  parameter::set(spars,"PHIMIN_DEG",PHIMIN_DEG_set);
  parameter::set(spars,"PHIMAX_DEG",PHIMAX_DEG_set);
  parameter::set(spars,"YMIN",YMIN_set);
  parameter::set(spars,"YMAX",YMAX_set);
}

void CSourceCalc_OSCAR::SetIDs(int IDa_set,int IDb_set){
  int IDa,IDb;
  parameter::set(spars,"IDa",IDa_set);
  parameter::set(spars,"IDb",IDb_set);
  IDa=IDa_set;
  IDb=IDb_set;
  if(abs(IDa)==211) parameter::set(spars,"Ma",139.58);
  else if(abs(IDa)==2212) parameter::set(spars,"Ma",938.28);
  else{
    printf("%d for IDa not recognized\n",IDa_set);
    exit(1);
  }

  if(abs(IDb)==211) parameter::set(spars,"Mb",139.58);
  else if(abs(IDb)==2212) parameter::set(spars,"Mb",938.28);
  else{
    printf("%d for IDb not recognized\n",IDb_set);
    exit(1);
  }
	  
  if(IDa==IDb) parameter::set(spars,"AEQUALB",1);

}

void CSourceCalc_OSCAR::CalcS(CCHArray *A){
  int alpha,nrmax;
  double cphi,sphi,ma,mb,phi,Pt,tau,volume;
  double rcm[4],**ra,**rb;
  int ir,ia,ib,na,nb,nbmax,IDa,IDb;
  double delr,x,y,z,xbar,ybar,zbar,x2bar,y2bar,z2bar,ex,ey,ez,snorm,r;
  bool AEQUALB;
  const double PI=4.0*atan(1.0);
  int NMAX=parameter::getI(spars,"NMAX",20000);
  nrmax=A->GetNRADIAL();
  delr=A->GetRADSTEP();
  
  printf("delr=%g, nrmax=%d\n",delr,nrmax);
  IDa=parameter::getI(spars,"IDa",211);
  IDb=parameter::getI(spars,"IDb",211);
  AEQUALB=parameter::getI(spars,"AEQUALB",0);
  if(IDa==IDb){
    AEQUALB=1;
    parameter::set(spars,"AEQUALB",1);
  }

  ra=new double *[NMAX];
  for(ia=0;ia<NMAX;ia++) ra[ia]=new double[4];
  if(AEQUALB){
    rb=ra;
  }
  else{
    rb=new double *[NMAX];
    for(ib=0;ib<NMAX;ib++) rb[ib]=new double[4];
  }
  ReadR(ra,rb,na,nb);

  xbar=ybar=zbar=x2bar=y2bar=z2bar=0.0;
  rcm[0]=0.0;
  for(ia=0;ia<na;ia++){
    if(AEQUALB) nbmax=ia-1;
    else nbmax=nb;
    for(ib=0;ib<nbmax;ib++){
      rcm[1]=ra[ia][1]-rb[ib][1];
      rcm[2]=ra[ia][2]-rb[ib][2];
      rcm[3]=ra[ia][3]-rb[ib][3];
      r=sqrt(rcm[1]*rcm[1]+rcm[2]*rcm[2]+rcm[3]*rcm[3]);
      x=rcm[1];
      y=rcm[2];
      z=rcm[3];
      xbar+=x;
      x2bar+=x*x;
      ybar+=y;
      y2bar+=y*y;
      zbar+=z;
      z2bar+=z*z;

      ir=int(floor(r/delr));
      if(ir<nrmax){
	ex=x/r; ey=y/r; ez=z/r;
	A->IncrementAExpArrayFromE(ex,ey,ez,1.0,ir);	
      }
    }
    if(10*(ia+1)%(10*int(na/10))==0)
      printf("finished %g percent\n",100*double(ia+1)/(10*int(na/10)));
  }
  A->FillRemainderX();

  if(AEQUALB) snorm=2.0/double(na*(na-1));
  else snorm=1.0/double(na*nb);
  xbar*=snorm;
  x2bar*=snorm;
  ybar*=snorm;
  y2bar*=snorm;
  zbar*=snorm;
  z2bar*=snorm;
  x2bar=x2bar-xbar*xbar;
  y2bar=y2bar-ybar*ybar;
  z2bar=z2bar-zbar*zbar;

  printf("xbar=%g, ybar=%g, zbar=%g\n",xbar,ybar,zbar);
  printf("Effective radii: Rout=%g, Rside=%g, Rlong=%g\n",
	 sqrt(0.5*x2bar),sqrt(0.5*y2bar),sqrt(0.5*z2bar));

  for(ir=0;ir<nrmax;ir++){
    volume=(4.0*PI/3)*(pow((ir+1)*delr,3)-pow(double(ir)*delr,3));
    A->ScaleArray(snorm/volume,ir);
  }

  for(ia=0;ia<NMAX;ia++){
    delete [] ra[ia];
  }
  delete [] ra;
  if(!AEQUALB){
    for(ib=0;ib<NMAX;ib++){
      delete [] rb[ib];
    }
    delete [] rb;
  }

}


void CSourceCalc_OSCAR::ReadR(double **ra,double **rb,int &na,int &nb){
  string OSCARfilename;
  ifstream oscarfile;
  double r[4],p[4];
  double Pt,DELPT,pt,pta,ptb,Ma,Mb,vperp;
  double phimin,phimax,phi,gamma,gammav;
  double YMIN,YMAX,y,sinhy,coshy,rout,rlong,rside,t,tau,mass;
  bool AEQUALB;
  int IDa,IDb,ident,idummy,i,ievent,nevents,NMAX,alpha;
  int npart,npartmax,ipart,identbonus;
  char dummy[160];
  const double PI=4.0*atan(1.0);
  const double TAUCOMPARE=12.0;
  na=nb=0;

  OSCARfilename=parameter::getS(spars,"OSCARfilename","UNDEFINED");
  nevents=parameter::getI(spars,"NEVENTSMAX",100);
  NMAX=parameter::getI(spars,"NMAX",10000);
  Pt=parameter::getD(spars,"Pt",600);
  DELPT=parameter::getD(spars,"DELPT",20);
  phimin=parameter::getD(spars,"PHIMIN_DEG",0);
  phimin=phimin*2.0*PI/360.0;
  phimax=parameter::getD(spars,"PHIMAX_DEG",360);
  phimax=phimax*2.0*PI/360.0;
  YMIN=parameter::getD(spars,"YMIN",-3.0);
  YMAX=parameter::getD(spars,"YMAX",3.0);
  Ma=parameter::getD(spars,"Ma",139.58);
  Mb=parameter::getD(spars,"Mb",139.58);
  pta=Pt*Ma/(Ma+Mb);
  ptb=Pt*Mb/(Ma+Mb);
  IDa=parameter::getI(spars,"IDa",211);
  IDb=parameter::getI(spars,"IDb",211);
  AEQUALB=parameter::getB(spars,"AEQUALB",bool(0));
  gammav=Pt/(sqrt(Ma*Ma+pta*pta)+sqrt(Mb*Mb+ptb*ptb));
  gammav=gammav/sqrt(1.0-gammav*gammav);
  gamma=sqrt(1.0+gammav*gammav);
  printf("Pt=%g, gammav=%g, gamma=%g\n",Pt,gammav,gamma);

  printf("Opening %s\n",OSCARfilename.c_str());
  oscarfile.open(OSCARfilename.c_str());

  for(i=0;i<5;i++){
    oscarfile.getline(dummy,80);
    printf("dummy line=%s\n",dummy);
  }


  //for(ievent=1;ievent<=nevents;ievent++){
  ievent=0;
  do{
    ievent+=1;
    oscarfile >> npartmax >> idummy;
    for(npart=0;npart<npartmax;npart++){
      oscarfile >> ipart >> ident >> identbonus >> p[1] >> p[2] >> p[3] 
		>> p[0] >> mass >> r[1] >> r[2] >> r[3] >> r[0];
      for(alpha=0;alpha<4;alpha++) p[alpha]*=1000.0;
      /* SET KINEMATIC VARIABLES */
      if(ident==IDa || ident==IDb){
	pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	if( (ident==IDa && fabs(pt-pta)<DELPT) ||
	    (ident==IDb && fabs(pt-ptb)<DELPT)){
	  phi=atan2(p[1],p[2]);
	  if(phi<0) phi+=2.0*PI;
	  if(phi>phimin && phi<phimax){
	    y=atanh(p[3]/p[0]);
	    if(y>YMIN && y<YMAX){
	  
	      if(ident==IDa) mass=Ma;
	      else mass=Mb;
	      p[0]=sqrt(pt*pt+p[3]*p[3]+mass*mass);
	      
	      rout=(p[1]*r[1]+p[2]*r[2])/pt;
	      rside=(p[1]*r[2]-p[2]*r[1])/pt;
	      sinhy=sinh(y);
	      coshy=cosh(y);
	      rlong=coshy*r[3]-sinhy*r[0];
	      tau=coshy*r[0]-sinhy*r[3];
	      vperp=pt/sqrt(mass*mass+pt*pt);
	      rout=rout-vperp*(tau-TAUCOMPARE);
	      tau=TAUCOMPARE;


	      //printf("r=(%g,%g,%g,%g)\n",tau,rout,rside,rlong);
	      if(ident==IDa || (ident==IDb && AEQUALB)){
		ra[na][0]=gamma*tau-gammav*rout;
		ra[na][1]=gamma*rout-gammav*tau;
		ra[na][2]=rside;
		ra[na][3]=rlong;
		na+=1;
		if(na==NMAX){
		  printf("Too many type-a particles, increase NMAX\n");
		  exit(1);
		}
	      }
	      else if(!AEQUALB){
		rb[nb][0]=gamma*tau-gammav*rout;
		rb[nb][1]=gamma*rout-gammav*tau;
		rb[nb][2]=rside;
		rb[nb][3]=rlong;
		nb+=1;
		if(nb==NMAX){
		  printf("Too many type-b particles, increase NMAX\n");
		  exit(1);
		}
		
	      }
	      
	    }
	  }  
	}
      }
    }
    oscarfile.getline(dummy,80);
    for(i=0;i<3;i++){
      if(!oscarfile.eof())
	oscarfile.getline(dummy,80);
      //printf("dummy line=%s\n",dummy);
    }
  } while(ievent<nevents&& !oscarfile.eof());
  oscarfile.close();
  if(AEQUALB) nb=na;
  printf("OSCAR file read: %d events, na=%d, nb=%d\n",ievent,na,nb);
}

#endif
