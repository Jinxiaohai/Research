void onepart_filter(int whichpart,double *p,int *ifilter){
  /* This should NOT depend on phi!
     If identical particles, this is only called with whichpart==1,
     If not identical, this gets called with whichpart==1 for the first
     particle and whichpart==2 for the second. */
  static double ptmax[20]={440.0,440.0,400.0,400.0,360.0,320.0,320.0,280.0,280.0,240.0,240.0,200.0,200.0,200.0,160.0,160.0,120.0,120.0,120.0,80.0};
  static double ptmin[20]={320.0,280.0,240.0,240.0,200.0,160.0,160.0,120.0,120.0,80.0,80.0,80.0,80.0,40.0,40.0,40.0,0.0,0.0,0.0,0.0};
  double rapidity,pt;
  int iy;
  *ifilter=1;
  rapidity=0.5*log((p[0]+p[3])/(p[0]-p[3]));
  rapidity=rapidity+2.92;
  if(rapidity<3.10 || rapidity>4.1){
    *ifilter=0;
    goto END_FILTER1 ;
  }
  iy=(int)floor(20.0*(rapidity-3.10));
  pt=sqrt(p[2]*p[2]+p[1]*p[1]);
  if(pt>ptmax[iy]||pt<ptmin[iy]){
    *ifilter=0;
    goto END_FILTER1 ;
  }
 END_FILTER1:;
}

/* ********************************************* */

/* These data were copied from John Sullivan's page:
   http://p2hp2.lanl.gov/people/sullivan/notes/na44/acceptance/
   ftn10_408_pbpb.dat.htm, May 12, 1998 */
void twopart_filter(double *p1,double *p2,int *ifilter){
  static double pycut[20][12]=
  {{5.53, 3.42, 0.06, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,},
   { 2.39,17.73, 7.64, 0.09, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,},
   { 0.00,10.48,18.81, 7.31, 0.04, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,},
   { 0.00, 0.78,17.37,17.20, 3.94, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,},
   { 0.00, 0.00, 3.82,19.53,14.68, 1.48, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,},
   { 0.00, 0.00, 0.00,13.13,17.65, 8.89, 0.11, 0.00, 0.00, 0.00, 0.00, 0.00,},
   { 0.00, 0.00, 0.00, 2.14,16.76,14.62, 3.40, 0.00, 0.00, 0.00, 0.00, 0.00,},
   { 0.00, 0.00, 0.00, 0.00, 9.92,16.40, 9.09, 0.13, 0.00, 0.00, 0.00, 0.00,},
   { 0.00, 0.00, 0.00, 0.00, 1.62,16.99,13.96, 3.07, 0.00, 0.00, 0.00, 0.00,},
   { 0.00, 0.00, 0.00, 0.00, 0.00,10.47,15.51, 8.27, 0.03, 0.00, 0.00, 0.00,},
   { 0.00, 0.00, 0.00, 0.00, 0.00, 3.13,17.09,11.73, 1.09, 0.00, 0.00, 0.00,},
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,13.06,14.21, 3.89, 0.00, 0.00, 0.00,},
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 6.30,15.59, 7.61, 0.00, 0.00, 0.00,},
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.39,14.78,12.55, 0.78, 0.00, 0.00,},
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 9.98,14.12, 3.98, 0.00, 0.00,},
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 3.02,15.51, 8.17, 0.00, 0.00,},
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,12.70,13.73, 0.85, 0.00,},
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 6.90,15.51, 7.51, 0.07,},
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.68,14.66,14.76, 3.51,},
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 3.93, 9.89, 0.72}};
  double pt44,rapidity;
  int iy,ipt;
  *ifilter=1;
  /* Check for first particle */
  rapidity=2.92+0.5*log((p1[0]+p1[3])/(p1[0]-p1[3]));
  iy=(int)floor(20.0*(rapidity-3.10));
  pt44=sqrt(p1[1]*p1[1]+p1[2]*p1[2]);
  if(p1[1]<0.0) pt44=-pt44;
  ipt=(int)floor((440.0-pt44)/40.0);
  if(ipt<0 || ipt>11){
    *ifilter=0;
    goto END_FILTER2;
  }
  if(fabs(p1[2])>pycut[iy][ipt]){
    *ifilter=0;
    goto END_FILTER2;
  }
  /* Check for second particle */
  rapidity=2.92+0.5*log((p2[0]+p2[3])/(p2[0]-p2[3]));
  iy=(int)floor(20.0*(rapidity-3.10));
  pt44=sqrt(p2[1]*p2[1]+p2[2]*p2[2]);
  if(p2[1]<0.0) pt44=-pt44;
  ipt=(int)floor((440.0-pt44)/40.0);
  if(ipt<0 || ipt>11){
    *ifilter=0;
    goto END_FILTER2;
  }
  if(fabs(p2[2])>pycut[iy][ipt]){
    *ifilter=0;
    goto END_FILTER2;
  }
 END_FILTER2:;
}
