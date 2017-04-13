        SUBROUTINE DECAY(IRUN,I,NNN,ISEED,wid,nt)
        PARAMETER (MAXSTR=150001,MAXR=1,
     1  AMN=0.939457,ETAM=0.5475,AMP=0.93828,AP1=0.13496,
     2  AP2=0.13957,AM0=1.232,PI=3.1415926)
        COMMON /AA/ R(3,MAXSTR)
        COMMON /BB/ P(3,MAXSTR)
        COMMON /CC/ E(MAXSTR)
        COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
        COMMON   /RUN/NUM
        COMMON   /PA/RPION(3,MAXSTR,MAXR)
        COMMON   /PB/PPION(3,MAXSTR,MAXR)
        COMMON   /PC/EPION(MAXSTR,MAXR)
        COMMON   /PD/LPION(MAXSTR,MAXR)
        COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &       IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
      COMMON/RNDF77/NSEED
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
      SAVE   
        lbanti=LB(I)
        DM=E(I)
        IF(iabs(LB(I)).EQ.11)THEN
           X3=RANART(NSEED)
           IF(X3.GT.(1./3.))THEN
              LB(I)=2
              NLAB=2
              LPION(NNN,IRUN)=5
              EPION(NNN,IRUN)=AP2
           ELSE
              LB(I)=1
              NLAB=1
              LPION(NNN,IRUN)=4
              EPION(NNN,IRUN)=AP1
           ENDIF
        ELSEIF(iabs(LB(I)).EQ.10)THEN
           X4=RANART(NSEED)
           IF(X4.GT.(1./3.))THEN
              LB(I)=1
              NLAB=1
              LPION(NNN,IRUN)=3
              EPION(NNN,IRUN)=AP2
           ELSE
              LB(I)=2
              NALB=2
              LPION(NNN,IRUN)=4
              EPION(NNN,IRUN)=AP1
           ENDIF
        ELSEIF(iabs(LB(I)).EQ.12)THEN
           CTRL=0.65
           IF(DM.lE.1.49)ctrl=-1.
           X5=RANART(NSEED)
           IF(X5.GE.ctrl)THEN
              X6=RANART(NSEED)
              IF(X6.GT.(1./3.))THEN
                 LB(I)=1
                 NLAB=1
                 LPION(NNN,IRUN)=3
                 EPION(NNN,IRUN)=AP2
              ELSE
                 LB(I)=2
                 NALB=2
                 LPION(NNN,IRUN)=4
                 EPION(NNN,IRUN)=AP1
              ENDIF
           ELSE
              LB(I)=2
              NLAB=2
              LPION(NNN,IRUN)=0
              EPION(NNN,IRUN)=ETAM
           ENDIF
        ELSEIF(iabs(LB(I)).EQ.13)THEN
           CTRL=0.65
           IF(DM.lE.1.49)ctrl=-1.
           X5=RANART(NSEED)
           IF(X5.GE.ctrl)THEN
              X8=RANART(NSEED)
              IF(X8.GT.(1./3.))THEN
                 LB(I)=2
                 NLAB=2
                 LPION(NNN,IRUN)=5
                 EPION(NNN,IRUN)=AP2
              ELSE
                 LB(I)=1
                 NLAB=1
                 LPION(NNN,IRUN)=4
                 EPION(NNN,IRUN)=AP1
              ENDIF
           ELSE
              LB(I)=1
              NLAB=1
              LPION(NNN,IRUN)=0
              EPION(NNN,IRUN)=ETAM
           ENDIF
        ENDIF
        CALL DKINE(IRUN,I,NNN,NLAB,ISEED,wid,nt)
        if(lbanti.lt.0) then
           lbi=LB(I)
           if(lbi.eq.1.or.lbi.eq.2) then
              lbi=-lbi
           elseif(lbi.eq.3) then
              lbi=5
           elseif(lbi.eq.5) then
              lbi=3
           endif
           LB(I)=lbi
           lbi=LPION(NNN,IRUN)
           if(lbi.eq.3) then
              lbi=5
           elseif(lbi.eq.5) then
              lbi=3
           elseif(lbi.eq.1.or.lbi.eq.2) then
              lbi=-lbi
           endif
           LPION(NNN,IRUN)=lbi
        endif
        if(nt.eq.ntmax) then
           lbm=LPION(NNN,IRUN)
           if(lbm.eq.0.or.lbm.eq.25
     1          .or.lbm.eq.26.or.lbm.eq.27) then
              lbsave=lbm
              xmsave=EPION(NNN,IRUN)
              pxsave=PPION(1,NNN,IRUN)
              pysave=PPION(2,NNN,IRUN)
              pzsave=PPION(3,NNN,IRUN)
              dpsave=dppion(NNN,IRUN)
              LPION(NNN,IRUN)=LB(I)
              EPION(NNN,IRUN)=E(I)
              PPION(1,NNN,IRUN)=P(1,I)
              PPION(2,NNN,IRUN)=P(2,I)
              PPION(3,NNN,IRUN)=P(3,I)
              dppion(NNN,IRUN)=dpertp(I)
              LB(I)=lbsave
              E(I)=xmsave
              P(1,I)=pxsave
              P(2,I)=pysave
              P(3,I)=pzsave
              dpertp(I)=dpsave
           endif
        endif
       RETURN
       END
