      SUBROUTINE ro2et2(i1,i2,lbb1,lbb2,ei1,ei2,iblock,iseed)
      PARAMETER (MAXSTR=150001)
      parameter (ETAM=0.5475,arho=0.77)
      COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
      COMMON/RNDF77/NSEED
      SAVE   
      if(lb(i1).ge.25.and.lb(i1).le.27.and.
     1     lb(i2).ge.25.and.lb(i2).le.27) then
         iblock=1895
         ei1=etam
         ei2=etam
         lbb1=0
         lbb2=0
      elseif(lb(i1).eq.0.and.lb(i2).eq.0) then
         iblock=1896
         lbb1=25+int(3*RANART(NSEED))
         lbb2=25+int(3*RANART(NSEED))
         ei1=arho
         ei2=arho
      endif
      return
      END
