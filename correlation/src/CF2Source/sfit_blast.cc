#ifndef __INCLUDE_SFIT_BLAST_CC__
#define __INCLUDE_SFIT_BLAST_CC__

CCF2SFit_Blast::CCF2SFit_Blast(CSourceCalc *scset,
			       C3DArray *cexpset,
			       C3DArray *cerrorset,
			       C3DArray *ctheory3Dset,
			       CCHArray *ctheoryset,
			       CCHArray *sourceset,
			       CKernel *kernelset){
  int i,j;
  ndim=3;

  // npars is different for different subclasses
  npars=4;
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

  par[0]->Set("lambda",parameter::getD(sourcecalc->spars,"lambda",1),
	     0.05,0.0,1.5);
  par[1]->Set("R",parameter::getD(sourcecalc->spars,"R",13),
	     0.5,5.0,20.0);
  par[2]->Set("Tau",parameter::getD(sourcecalc->spars,"Tau",12),
	     0.5,5.0,20.0);
  par[3]->Set("DelTau",parameter::getD(sourcecalc->spars,"DelTau",5),
	     0.5,0.0,20.0);


  for(i=0;i<npars;i++){
    for(j=0;j<npars;j++){
      StepMatrix[i][j]=0.0;
      ErrorMatrix[i][j]=0.0;
      if(i==j){
	StepMatrix[i][j]=par[i]->error;
	ErrorMatrix[i][j]=par[i]->error*par[i]->error;
      }
    }
  }
}
#endif

