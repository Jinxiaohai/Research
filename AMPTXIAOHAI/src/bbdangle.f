      subroutine bbdangle(pxd,pyd,pzd,nt,ipert1,ianti,idloop,pfinal,
     1 dprob1,lbm)
      PARAMETER (PI=3.1415926)
      common/leadng/lb1,px1,py1,pz1,em1,e1,xfnl,yfnl,zfnl,tfnl,
     1     px1n,py1n,pz1n,dp1n
      common /dpi/em2,lb2
      COMMON/RNDF77/NSEED
      common /para8/ idpert,npertd,idxsec
      COMMON /AREVT/ IAEVT, IARUN, MISS
      SAVE   
      C1=1.0-2.0*RANART(NSEED)
      T1=2.0*PI*RANART(NSEED)
      S1=SQRT(1.0-C1**2)
      CT1=COS(T1)
      ST1=SIN(T1)
      PZd=pfinal*C1
      PXd=pfinal*S1*CT1 
      PYd=pfinal*S1*ST1
      if(idpert.eq.1.and.npertd.ge.1) then
         dprob=dprob1
      elseif(idpert.eq.2.and.npertd.ge.1) then
         dprob=1./float(npertd)
      endif
      if(ianti.eq.0) then
         if(idpert.eq.0.or.(idpert.eq.1.and.ipert1.eq.0).or.
     1        (idpert.eq.2.and.idloop.eq.(npertd+1))) then
            write (91,*) lb1,' *',lb2,' ->d+',lbm,' (regular d prodn) 
     1 @evt#',iaevt,' @nt=',nt
         elseif((idpert.eq.1.or.idpert.eq.2).and.idloop.eq.npertd) then
            write (91,*) lb1,' *',lb2,' ->d+',lbm,' (pert d prodn) 
     1 @evt#',iaevt,' @nt=',nt,' @prob=',dprob
         endif
      else
         if(idpert.eq.0.or.(idpert.eq.1.and.ipert1.eq.0).or.
     1        (idpert.eq.2.and.idloop.eq.(npertd+1))) then
            write (91,*) lb1,' *',lb2,' ->d+',lbm,' (regular dbar prodn) 
     1 @evt#',iaevt,' @nt=',nt
         elseif((idpert.eq.1.or.idpert.eq.2).and.idloop.eq.npertd) then
            write (91,*) lb1,' *',lb2,' ->d+',lbm,' (pert dbar prodn) 
     1 @evt#',iaevt,' @nt=',nt,' @prob=',dprob
         endif
      endif
      return
      end
