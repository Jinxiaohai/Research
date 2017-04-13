      SUBROUTINE rpiret(i1,i2,lbb1,lbb2,ei1,ei2,iblock,iseed)
      PARAMETER (MAXSTR=150001)
      PARAMETER (AP1=0.13496,AP2=0.13957,ETAM=0.5475,arho=0.77)
      COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
      COMMON/RNDF77/NSEED
      SAVE   
      if((lb(i1).ge.25.and.lb(i1).le.27
     1     .and.lb(i2).ge.3.and.lb(i2).le.5).or.
     2     (lb(i1).ge.3.and.lb(i1).le.5
     3     .and.lb(i2).ge.25.and.lb(i2).le.27)) then
         iblock=1880
         ei1=arho
         ei2=etam
         lbb1=25+int(3*RANART(NSEED))
         lbb2=0
      elseif((lb(i1).ge.25.and.lb(i1).le.27.and.lb(i2).eq.0).or.
     1        (lb(i2).ge.25.and.lb(i2).le.27.and.lb(i1).eq.0)) then
         iblock=1881
         lbb1=25+int(3*RANART(NSEED))
         lbb2=3+int(3*RANART(NSEED))
         ei1=arho
         ei2=ap2
         if(lbb2.eq.4) ei2=ap1
      endif
      return
      END
