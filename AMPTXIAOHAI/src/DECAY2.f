        SUBROUTINE DECAY2(IRUN,I,NNN,ISEED,wid,nt)
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
      COMMON/RNDF77/NSEED
      SAVE   
        lbanti=LB(I)
        DM=E(I)
        IF(iabs(LB(I)).EQ.11)THEN
           X3=RANART(NSEED)
           IF(X3.LT.(1./3))THEN
              LB(I)=2
              NLAB=2
              LPION(NNN,IRUN)=5
              EPION(NNN,IRUN)=AP2
              LPION(NNN+1,IRUN)=4
              EPION(NNN+1,IRUN)=AP1
           ELSEIF(X3.LT.2./3.AND.X3.GT.1./3.)THEN
              LB(I)=1
              NLAB=1
              LPION(NNN,IRUN)=5
              EPION(NNN,IRUN)=AP2
              LPION(NNN+1,IRUN)=3
              EPION(NNN+1,IRUN)=AP2
           ELSE
              LB(I)=1
              NLAB=1
              LPION(NNN,IRUN)=4
              EPION(NNN,IRUN)=AP1
              LPION(NNN+1,IRUN)=4
              EPION(NNN+1,IRUN)=AP1
           ENDIF
        ELSEIF(iabs(LB(I)).EQ.10)THEN
           X3=RANART(NSEED)
           IF(X3.LT.(1./3))THEN
              LB(I)=2
              NLAB=2
              LPION(NNN,IRUN)=4
              EPION(NNN,IRUN)=AP1
              LPION(NNN+1,IRUN)=4
              EPION(NNN+1,IRUN)=AP1
           ELSEIF(X3.LT.2./3.AND.X3.GT.1./3.)THEN
              LB(I)=1
              NLAB=1
              LPION(NNN,IRUN)=3
              EPION(NNN,IRUN)=AP2
              LPION(NNN+1,IRUN)=4
              EPION(NNN+1,IRUN)=AP1
           ELSE
              LB(I)=2
              NLAB=2
              LPION(NNN,IRUN)=5
              EPION(NNN,IRUN)=AP2
              LPION(NNN+1,IRUN)=3
              EPION(NNN+1,IRUN)=AP2
           ENDIF
        ENDIF
        CALL DKINE2(IRUN,I,NNN,NLAB,ISEED,wid,nt)
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
           lbi=LPION(NNN+1,IRUN)
           if(lbi.eq.3) then
              lbi=5
           elseif(lbi.eq.5) then
              lbi=3
           elseif(lbi.eq.1.or.lbi.eq.2) then
              lbi=-lbi
           endif
           LPION(NNN+1,IRUN)=lbi
        endif
       RETURN
       END
