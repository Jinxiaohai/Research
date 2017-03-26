#ifndef __INCLUDE_SFIT_METROPOLIS_3DGAUSSIAN_CC__
#define __INCLUDE_SFIT_METROPOLIS_3DGAUSSIAN_CC__

CCF2S_Metropolis_3DGaussian::CCF2S_Metropolis_3DGaussian(CSourceCalc *scset,
							 C3DArray *cexpset,
							 C3DArray *cerrorset,
							 C3DArray *ctheory3Dset,
							 CCHArray *ctheoryset,
							 CCHArray *sourceset,
							 CKernel *kernelset){
  int i,j;
  ndim=3;

  // npars is different for different subclasses
  npars=5;
  //
  sourcecalc=scset;
  cexp3D=cexpset;
  cerror3D=cerrorset;
  ctheory=ctheoryset;
  ctheory3D=ctheory3Dset;
  source=sourceset;
  kernel=kernelset;
  Init();

  // initialization of pars is also unique to given subclass
  par[0].Set("lambda",parameter::getD(sourcecalc->spars,"lambda",1),
	     0.02,0.0,1.5);
  par[1].Set("Rx",parameter::getD(sourcecalc->spars,"Rx",5),
	     0.2,1.0,10.0);
  par[2].Set("Ry",parameter::getD(sourcecalc->spars,"Ry",5),
	     0.2,1.0,10.0);
  par[3].Set("Rz",parameter::getD(sourcecalc->spars,"Rz",5),
	     0.2,1.0,10.0);
  par[4].Set("Xoff",parameter::getD(sourcecalc->spars,"Xoff",0),
	     0.2,-10.0,10.0);
  for(i=0;i<npars;i++){
    for(j=0;j<npars;j++){
      Usigma[i][j]=0.0;
      if(i==j) Usigma[i][j]=par[i].error;
    }
  }

}
#endif

