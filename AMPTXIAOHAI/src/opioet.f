      SUBROUTINE opioet(i1,i2,lbb1,lbb2,ei1,ei2,iblock,iseed)
      PARAMETER (MAXSTR=150001)
      PARAMETER (AP1=0.13496,AP2=0.13957,ETAM=0.5475,aomega=0.782)
      COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
      COMMON/RNDF77/NSEED
      SAVE   
      if((lb(i1).ge.3.and.lb(i1).le.5.and.lb(i2).eq.28).or.
     1     (lb(i2).ge.3.and.lb(i2).le.5.and.lb(i1).eq.28)) then
         iblock=1890
         ei1=aomega
         ei2=etam
         lbb1=28
         lbb2=0
      elseif((lb(i1).eq.28.and.lb(i2).eq.0).or.
     1        (lb(i1).eq.0.and.lb(i2).eq.28)) then
         iblock=1891
         lbb1=28
         lbb2=3+int(3*RANART(NSEED))
         ei1=aomega
         ei2=ap2
         if(lbb2.eq.4) ei2=ap1
      endif
      return
      END
