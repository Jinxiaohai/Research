      SUBROUTINE CRPP(PX,PY,PZ,SRT,I1,I2,IBLOCK,
     &ppel,ppin,spprho,ipp)
      PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1     AMP=0.93828,AP1=0.13496,
     2 AP2=0.13957,AM0=1.232,PI=3.1415926,CUTOFF=1.8966,AVMASS=0.9383)
      PARAMETER      (AKA=0.498,aks=0.895)
      parameter     (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
      COMMON /AA/ R(3,MAXSTR)
      COMMON /BB/ P(3,MAXSTR)
      COMMON /CC/ E(MAXSTR)
      COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
      COMMON/RNDF77/NSEED
      SAVE   
      lb1i=lb(i1)
      lb2i=lb(i2)
       PX0=PX
       PY0=PY
       PZ0=PZ
        iblock=1
        if(srt.gt.(2*aka).and.(ppin/(ppin+ppel)).gt.RANART(NSEED)) then
           ranpi=RANART(NSEED)
           if((pprr/ppin).ge.ranpi) then
              call pi2ro2(i1,i2,lbb1,lbb2,ei1,ei2,iblock,iseed)
           elseif((pprr+ppee)/ppin.ge.ranpi) then
              call pi2et2(i1,i2,lbb1,lbb2,ei1,ei2,iblock,iseed)
           elseif(((pprr+ppee+pppe)/ppin).ge.ranpi) then
              call pi3eta(i1,i2,lbb1,lbb2,ei1,ei2,iblock,iseed)
           elseif(((pprr+ppee+pppe+rpre)/ppin).ge.ranpi) then
              call rpiret(i1,i2,lbb1,lbb2,ei1,ei2,iblock,iseed)
           elseif(((pprr+ppee+pppe+rpre+xopoe)/ppin).ge.ranpi) then
              call opioet(i1,i2,lbb1,lbb2,ei1,ei2,iblock,iseed)
           elseif(((pprr+ppee+pppe+rpre+xopoe+rree)
     1             /ppin).ge.ranpi) then
              call ro2et2(i1,i2,lbb1,lbb2,ei1,ei2,iblock,iseed)
           elseif(((pprr+ppee+pppe+rpre+xopoe+rree+ppinnb)/ppin)
     1             .ge.ranpi) then
              call bbarfs(lbb1,lbb2,ei1,ei2,iblock,iseed)
           else
              iblock=66
              ei1=aka
              ei2=aka
              lbb1=21
              lbb2=23
              lb1=lb(i1)
              lb2=lb(i2)
        if( ( (lb1.eq.0.or.(lb1.ge.3.and.lb1.le.5))
     1       .and.(lb2.ge.25.and.lb2.le.28))
     2       .or. ( (lb2.eq.0.or.(lb2.ge.3.and.lb2.le.5))
     3       .and.(lb1.ge.25.and.lb1.le.28))) then
           ei1=aks
           ei2=aka
           if(RANART(NSEED).ge.0.5) then
              iblock=366
              lbb1=30
              lbb2=21
           else
              iblock=367
              lbb1=-30
              lbb2=23
           endif
        endif
           endif
           e(i1)=ei1
           e(i2)=ei2
           lb(i1)=lbb1
           lb(i2)=lbb2
       else
         if ((lb(i1).lt.3.or.lb(i1).gt.5).and.
     &        (lb(i2).lt.3.or.lb(i2).gt.5)) return
        IBLOCK=6
       if(ipp.eq.1.or.ipp.eq.4.or.ipp.eq.6)go to 10
       if(spprho/ppel.gt.RANART(NSEED))go to 20
       endif
10      NTAG=0
        EM1=E(I1)
        EM2=E(I2)
          PR2   = (SRT**2 - EM1**2 - EM2**2)**2
     1                - 4.0 * (EM1*EM2)**2
          IF(PR2.LE.0.)PR2=1.e-09
          PR=SQRT(PR2)/(2.*SRT)
          C1   = 1.0 - 2.0 * RANART(NSEED)
          T1   = 2.0 * PI * RANART(NSEED)
          S1   = SQRT( 1.0 - C1**2 )
      CT1  = COS(T1)
      ST1  = SIN(T1)
      PZ   = PR * C1
      PX   = PR * S1*CT1 
      PY   = PR * S1*ST1
      CALL ROTATE(PX0,PY0,PZ0,PX,PY,PZ) 
      RETURN
20       continue
       iblock=666
       call rhores(i1,i2)
       if(ipp.eq.2)lb(i1)=27
       if(ipp.eq.3)lb(i1)=26
       if(ipp.eq.5)lb(i1)=25
       return       
      END
