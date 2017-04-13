      subroutine dmangle(pxn,pyn,pzn,nt,ianti,pfinal,lbm)
      PARAMETER (PI=3.1415926)
      common/leadng/lb1,px1,py1,pz1,em1,e1,xfnl,yfnl,zfnl,tfnl,
     1     px1n,py1n,pz1n,dp1n
      common /dpi/em2,lb2
      COMMON /AREVT/ IAEVT, IARUN, MISS
      COMMON/RNDF77/NSEED
      SAVE   
      C1=1.0-2.0*RANART(NSEED)
      T1=2.0*PI*RANART(NSEED)
      S1=SQRT(1.0-C1**2)
      CT1=COS(T1)
      ST1=SIN(T1)
      Pzn=pfinal*C1
      Pxn=pfinal*S1*CT1 
      Pyn=pfinal*S1*ST1
      if(ianti.eq.0) then
         write (91,*) ' d+',lbm,' ->BB (regular d destrn) @evt#',
     1        iaevt,' @nt=',nt,' lb1,2=',lb1,lb2
      else
         write (91,*) ' d+',lbm,' ->BB (regular dbar destrn) @evt#',
     1        iaevt,' @nt=',nt,' lb1,2=',lb1,lb2
      endif
      return
      end
