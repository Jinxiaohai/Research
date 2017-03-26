void momentum_smear(double *p1,double *p2,double *mom,double *qred){
#if ! defined IDENTICAL
  static double pdotk=MASS2*MASS2-MASS1*MASS1;
#endif
  double delpg[4],delp,pmag2;
  int alpha;
  /* smear the momenta by delp */
  
  delp=10.0;
  gauss(&delpg[0],&delpg[1]);
  gauss(&delpg[2],&delpg[3]);
  pmag2=0.0;
  for(alpha=1;alpha<4;alpha++){
    p1[alpha]=p1[alpha]+delp*delpg[alpha];
    pmag2=pmag2+(p1[alpha]*p1[alpha]);
  }
  p1[0]=sqrt(MASS1*MASS1+pmag2);

  //delp=10.0;
  gauss(&delpg[0],&delpg[1]);
  gauss(&delpg[2],&delpg[3]);
  pmag2=0.0;
  for(alpha=1;alpha<4;alpha++){
    p2[alpha]=p2[alpha]+delp*delpg[alpha];
    pmag2=pmag2+(p2[alpha]*p2[alpha]);
  }
  p2[0]=sqrt(MASS2*MASS2+pmag2);

#ifdef IDENTICAL
  *mom=-(p2[0]-p1[0])*(p2[0]-p1[0]);
  for(alpha=1;alpha<4;alpha++){
    *mom=*mom+(p2[alpha]-p1[alpha])*(p2[alpha]-p1[alpha]);
  }
  *mom=sqrt(*mom);
#else
  pmag2=(p1[0]+p2[0])*(p1[0]+p2[0]);
  *mom=-(p2[0]-p1[0])*(p2[0]-p1[0]);
  for(alpha=1;alpha<4;alpha++){
    pmag2=pmag2-(p1[alpha]+p2[alpha])*(p1[alpha]+p2[alpha]);
    *mom=*mom+pow(p2[alpha]-p1[alpha])*(p2[alpha]-p1[alpha]);
  }
  *mom=sqrt(*mom+pdotk*pdotk/(pmag2));
#endif
  *qred=0.5**mom;

}
