      SUBROUTINE CRND(IRUN,PX,PY,PZ,SRT,I1,I2,IBLOCK,
     &SIGNN,SIG,sigk,xsk1,xsk2,xsk3,xsk4,xsk5,NT,ipert1)
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,AKA=0.498,APHI=1.020,
     2  AP2=0.13957,AM0=1.232,PI=3.1415926,CUTOFF=1.8966,AVMASS=0.9383)
        parameter     (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
        parameter (xmd=1.8756,npdmax=10000)
        COMMON /AA/ R(3,MAXSTR)
        COMMON /BB/ P(3,MAXSTR)
        COMMON /CC/ E(MAXSTR)
        COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
        common /ff/f(-mx:mx,-my:my,-mz:mz,-mpx:mpx,-mpy:mpy,-mpz:mpzp)
        common /gg/ dx,dy,dz,dpx,dpy,dpz
        COMMON /INPUT/ NSTAR,NDIRCT,DIR
        COMMON /NN/NNN
        COMMON /BG/BETAX,BETAY,BETAZ,GAMMA
        COMMON   /RUN/NUM
        COMMON   /PA/RPION(3,MAXSTR,MAXR)
        COMMON   /PB/PPION(3,MAXSTR,MAXR)
        COMMON   /PC/EPION(MAXSTR,MAXR)
        COMMON   /PD/LPION(MAXSTR,MAXR)
        common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      common/leadng/lb1,px1,py1,pz1,em1,e1,xfnl,yfnl,zfnl,tfnl,
     1 px1n,py1n,pz1n,dp1n
      COMMON/RNDF77/NSEED
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
      common /dpi/em2,lb2
      common /para8/ idpert,npertd,idxsec
      dimension ppd(3,npdmax),lbpd(npdmax)
      SAVE   
       n12=0
       m12=0
        IBLOCK=0
        NTAG=0
        EM1=E(I1)
        EM2=E(I2)
        PR  = SQRT( PX**2 + PY**2 + PZ**2 )
        C2  = PZ / PR
        X1  = RANART(NSEED)
        ianti=0
        if(lb(i1).lt.0 .and. lb(i2).lt.0)ianti=1
      call sbbdm(srt,sdprod,ianti,lbm,xmm,pfinal)
      if(idpert.eq.1.and.ipert1.eq.1) then
         IF (SRT .LT. 2.012) RETURN
         if((iabs(lb(i1)).eq.1.or.iabs(lb(i1)).eq.2)
     1        .and.(iabs(lb(i2)).ge.6.and.iabs(lb(i2)).le.13)) then
            goto 108
         elseif((iabs(lb(i2)).eq.1.or.iabs(lb(i2)).eq.2)
     1           .and.(iabs(lb(i1)).ge.6.and.iabs(lb(i1)).le.13)) then
            goto 108
         else
            return
         endif
      endif
      IF (X1 .LE. SIGNN/SIG) THEN
        AS  = ( 3.65 * (SRT - 1.8766) )**6
        A   = 6.0 * AS / (1.0 + AS)
        TA  = -2.0 * PR**2
        X   = RANART(NSEED)
        T1  = sngl(DLOG(dble(1.-X)*DEXP(dble(A)*dble(TA))+dble(X)))/  A
        C1  = 1.0 - T1/TA
        T1  = 2.0 * PI * RANART(NSEED)
        IBLOCK=1
       GO TO 107
      ELSE
        IF (SRT .LT. 2.04) RETURN
        if(((iabs(LB(I1)).EQ.2.or.iabs(LB(I2)).EQ.2).AND.
     1       (LB(I1)*LB(I2)).EQ.20).or.(LB(I1)*LB(I2)).EQ.13) then
           IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
        ENDIF
        PRF=SQRT(0.25*SRT**2-AVMASS**2)
        IF(EM1.GT.1.)THEN
        DELTAM=EM1
        ELSE
        DELTAM=EM2
        ENDIF
        RENOM=DELTAM*PRF**2/DENOM(SRT,1.)/PR
        RENOMN=DELTAM*PRF**2/DENOM(SRT,2.)/PR
        RENOM1=DELTAM*PRF**2/DENOM(SRT,-1.)/PR
       if((iabs(lb(i1)).eq.2).and.(iabs(lb(i2)).eq.6)) renom=0.
       if((iabs(lb(i2)).eq.2).and.(iabs(lb(i1)).eq.6)) renom=0.
       if((iabs(lb(i1)).eq.1).and.(iabs(lb(i2)).eq.9)) renom=0.
       if((iabs(lb(i2)).eq.1).and.(iabs(lb(i1)).eq.9)) renom=0.
       Call M1535(iabs(lb(i1)),iabs(lb(i2)),srt,x1535)
        X1440=(3./4.)*SIGMA(SRT,2,0,1)
       if(((iabs(lb(i1)).eq.2).and.(iabs(lb(i2)).eq.6)).OR. 
     &         ((iabs(lb(i2)).eq.2).and.(iabs(lb(i1)).eq.6)).OR.
     &         ((iabs(lb(i1)).eq.1).and.(iabs(lb(i2)).eq.9)).OR.
     &         ((iabs(lb(i2)).eq.1).and.(iabs(lb(i1)).eq.9)))THEN
          IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
          IF((SIGK+SIGNN+sdprod)/SIG.GE.X1)GO TO 306
       ENDIF
        IF(LB(I1)*LB(I2).EQ.18.AND.
     &  (iabs(LB(I1)).EQ.2.OR.iabs(LB(I2)).EQ.2))then
        SIGND=SIGMA(SRT,1,1,0)+0.5*SIGMA(SRT,1,1,1)
        SIGDN=0.25*SIGND*RENOM
        IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
        IF(X1.GT.(SIGNN+SIGDN+X1440+X1535+SIGK+sdprod)/SIG)RETURN
       IF(SIGK/(SIGK+SIGDN+X1440+X1535).GT.RANART(NSEED))GO TO 306
       IF(RANART(NSEED).LT.SIGDN/(SIGDN+X1440+X1535))THEN
        M12=3
       GO TO 206
       ELSE
              IF(RANART(NSEED).LT.X1440/(X1440+X1535))THEN
              M12=37
              ELSE
                   return
              ENDIF
       GO TO 204
       ENDIF
        ENDIF
        IF(LB(I1)*LB(I2).EQ.6.AND.
     &   ((iabs(LB(I1)).EQ.1).OR.(iabs(LB(I2)).EQ.1)))then
        SIGND=SIGMA(SRT,1,1,0)+0.5*SIGMA(SRT,1,1,1)
        SIGDN=0.25*SIGND*RENOM
        IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
        IF (X1.GT.(SIGNN+SIGDN+X1440+X1535+SIGK+sdprod)/SIG)RETURN
       IF(SIGK/(SIGK+SIGDN+X1440+X1535).GT.RANART(NSEED))GO TO 306
       IF(RANART(NSEED).LT.SIGDN/(SIGDN+X1440+X1535))THEN
        M12=6
       GO TO 206
       ELSE
              IF(RANART(NSEED).LT.X1440/(X1440+X1535))THEN
              M12=47
              ELSE
                   return
              ENDIF
       GO TO 204
       ENDIF
        ENDIF
        IF(LB(I1)*LB(I2).EQ.8.AND.
     &   (iabs(LB(I1)).EQ.1.OR.iabs(LB(I2)).EQ.1))THEN
        SIGND=1.5*SIGMA(SRT,1,1,1)
        SIGDN=0.25*SIGND*RENOM
        IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
        IF(X1.GT.(SIGNN+SIGDN+x1440+x1535+SIGK+sdprod)/SIG)RETURN
       IF(SIGK/(SIGK+SIGDN+X1440+X1535).GT.RANART(NSEED))GO TO 306
       IF(RANART(NSEED).LT.SIGDN/(SIGDN+X1440+X1535))THEN
        M12=4
       GO TO 206
       ELSE
              IF(RANART(NSEED).LT.X1440/(X1440+X1535))THEN
              M12=39
              ELSE
              M12=40
              ENDIF
              GO TO 204
       ENDIF
        ENDIF
        IF(LB(I1)*LB(I2).EQ.14.AND.
     &   (iabs(LB(I1)).EQ.2.OR.iabs(LB(I2)).EQ.2))THEN
        SIGND=1.5*SIGMA(SRT,1,1,1)
        SIGDN=0.25*SIGND*RENOM
        IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
        IF(X1.GT.(SIGNN+SIGDN+x1440+x1535+SIGK+sdprod)/SIG)RETURN
       IF(SIGK/(SIGK+SIGDN+X1440+X1535).GT.RANART(NSEED))GO TO 306
       IF(RANART(NSEED).LT.SIGDN/(SIGDN+X1440+X1535))THEN
        M12=5
       GO TO 206
       ELSE
              IF(RANART(NSEED).LT.X1440/(X1440+X1535))THEN
              M12=48
              ELSE
              M12=49
              ENDIF
              GO TO 204
       ENDIF
        ENDIF
        IF(LB(I1)*LB(I2).EQ.16.AND.
     &   (iabs(LB(I1)).EQ.2.OR.iabs(LB(I2)).EQ.2))THEN
        SIGND=0.5*SIGMA(SRT,1,1,1)+0.25*SIGMA(SRT,1,1,0)
        SIGDN=0.5*SIGND*RENOM
        IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
        IF(X1.GT.(SIGNN+SIGDN+2.*x1440+2.*x1535+SIGK+sdprod)/SIG)RETURN
       IF(SIGK/(SIGK+SIGDN+2*X1440+2*X1535).GT.RANART(NSEED))GO TO 306
       IF(RANART(NSEED).LT.SIGDN/(SIGDN+2.*X1440+2.*X1535))THEN
        M12=1
       GO TO 206
       ELSE
              IF(RANART(NSEED).LT.X1440/(X1440+X1535))THEN
              M12=41
              IF(RANART(NSEED).LE.0.5)M12=43
              ELSE
              M12=42
              IF(RANART(NSEED).LE.0.5)M12=44
              ENDIF
              GO TO 204
       ENDIF
        ENDIF
        IF(LB(I1)*LB(I2).EQ.7)THEN
        SIGND=0.5*SIGMA(SRT,1,1,1)+0.25*SIGMA(SRT,1,1,0)
        SIGDN=0.5*SIGND*RENOM
        IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
        IF(X1.GT.(SIGNN+SIGDN+2.*x1440+2.*x1535+SIGK+sdprod)/SIG)RETURN
       IF(SIGK/(SIGK+SIGDN+2*X1440+2*X1535).GT.RANART(NSEED))GO TO 306
       IF(RANART(NSEED).LT.SIGDN/(SIGDN+2.*X1440+2.*X1535))THEN
        M12=2
       GO TO 206
       ELSE
              IF(RANART(NSEED).LT.X1440/(X1440+X1535))THEN
              M12=50
              IF(RANART(NSEED).LE.0.5)M12=51
              ELSE
              M12=52
              IF(RANART(NSEED).LE.0.5)M12=53
              ENDIF
              GO TO 204
       ENDIF
        ENDIF
        IF(LB(I1)*LB(I2).EQ.10.AND.
     &  (iabs(LB(I1)).EQ.1.OR.iabs(LB(I2)).EQ.1))then
        SIGND=(3./4.)*SIGMA(SRT,2,0,1)
        SIGDN=SIGND*RENOMN
        IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
        IF(X1.GT.(SIGNN+SIGDN+X1535+SIGK+sdprod)/SIG)RETURN
       IF(SIGK/(SIGK+SIGDN+X1535).GT.RANART(NSEED))GO TO 306
       IF(RANART(NSEED).LT.SIGDN/(SIGDN+X1535))THEN
        M12=7
        GO TO 206
       ELSE
       M12=54
       IF(RANART(NSEED).LE.0.5)M12=55
       ENDIF
       GO TO 204
        ENDIF
        IF(LB(I1)*LB(I2).EQ.22.AND.
     &   (iabs(LB(I1)).EQ.2.OR.iabs(LB(I2)).EQ.2))then
        SIGND=(3./4.)*SIGMA(SRT,2,0,1)
        SIGDN=SIGND*RENOMN
        IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
        IF(X1.GT.(SIGNN+SIGDN+X1535+SIGK+sdprod)/SIG)RETURN
       IF(SIGK/(SIGK+SIGDN+X1535).GT.RANART(NSEED))GO TO 306
       IF(RANART(NSEED).LT.SIGDN/(SIGDN+X1535))THEN
        M12=8
        GO TO 206
       ELSE
       M12=56
       IF(RANART(NSEED).LE.0.5)M12=57
       ENDIF
       GO TO 204
        ENDIF
        IF((iabs(LB(I1)).EQ.12).OR.(iabs(LB(I1)).EQ.13).OR.
     1  (iabs(LB(I2)).EQ.12).OR.(iabs(LB(I2)).EQ.13))THEN
        SIGND=X1535
        SIGDN=SIGND*RENOM1
        IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
        IF(X1.GT.(SIGNN+SIGDN+SIGK+sdprod)/SIG)RETURN
       IF(SIGK/(SIGK+SIGDN).GT.RANART(NSEED))GO TO 306
        IF(LB(I1)*LB(I2).EQ.24)M12=10
        IF(LB(I1)*LB(I2).EQ.12)M12=12
        IF(LB(I1)*LB(I2).EQ.26)M12=11
       IF(LB(I1)*LB(I2).EQ.13)M12=9
       GO TO 206
        ENDIF
204       CONTINUE
          DMAX = SRT - AVMASS-0.005
          DMIN = 1.078
          IF((M12.eq.37).or.(M12.eq.39).or.
     1    (M12.eQ.41).OR.(M12.eQ.43).OR.(M12.EQ.46).
     2     OR.(M12.EQ.48).OR.(M12.EQ.50).OR.(M12.EQ.51))then
          IF(DMAX.LT.1.44) THEN
          FM=FNS(DMAX,SRT,0.)
          ELSE
             xdmass=1.44
          FM=FNS(xdmass,SRT,1.)
          ENDIF
          IF(FM.EQ.0.)FM=1.E-09
          NTRY2=0
11        DM=RANART(NSEED)*(DMAX-DMIN)+DMIN
          NTRY2=NTRY2+1
          IF((RANART(NSEED).GT.FNS(DM,SRT,1.)/FM).AND.
     1    (NTRY2.LE.10)) GO TO 11
          if(dm.gt.2.14) goto 11
              GO TO 13
              ELSE
          IF(DMAX.LT.1.535) THEN
          FM=FD5(DMAX,SRT,0.)
          ELSE
             xdmass=1.535
          FM=FD5(xdmass,SRT,1.)
          ENDIF
          IF(FM.EQ.0.)FM=1.E-09
          NTRY1=0
12        DM = RANART(NSEED) * (DMAX-DMIN) + DMIN
          NTRY1=NTRY1+1
          IF((RANART(NSEED) .GT. FD5(DM,SRT,1.)/FM).AND.
     1    (NTRY1.LE.10)) GOTO 12
          if(dm.gt.1.84) goto 12
             ENDIF
13       CONTINUE
       PRF=0.
       PF2=((SRT**2-DM**2+AVMASS**2)/(2.*SRT))**2-AVMASS**2
       IF(PF2.GT.0.)PRF=SQRT(PF2)
          IF(M12.EQ.37)THEN
          IF(iabs(LB(I1)).EQ.9)THEN
          LB(I1)=1
          E(I1)=AMP
         LB(I2)=11
         E(I2)=DM
          ELSE
          LB(I2)=1
          E(I2)=AMP
         LB(I1)=11
         E(I1)=DM
          ENDIF
         GO TO 207
          ENDIF
          IF(M12.EQ.38)THEN
          IF(iabs(LB(I1)).EQ.9)THEN
          LB(I1)=1
          E(I1)=AMP
         LB(I2)=13
         E(I2)=DM
          ELSE
          LB(I2)=1
          E(I2)=AMP
         LB(I1)=13
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
          IF(M12.EQ.39)THEN
          IF(iabs(LB(I1)).EQ.8)THEN
          LB(I1)=1
          E(I1)=AMP
         LB(I2)=11
         E(I2)=DM
          ELSE
          LB(I2)=1
          E(I2)=AMP
         LB(I1)=11
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
          IF(M12.EQ.40)THEN
          IF(iabs(LB(I1)).EQ.8)THEN
          LB(I1)=1
          E(I1)=AMP
         LB(I2)=13
         E(I2)=DM
          ELSE
          LB(I2)=1
          E(I2)=AMP
         LB(I1)=13
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
          IF(M12.EQ.41)THEN
          IF(iabs(LB(I1)).EQ.8)THEN
          LB(I1)=2
          E(I1)=AMN
         LB(I2)=11
         E(I2)=DM
          ELSE
          LB(I2)=2
          E(I2)=AMN
         LB(I1)=11
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
          IF(M12.EQ.42)THEN
          IF(iabs(LB(I1)).EQ.8)THEN
          LB(I1)=2
          E(I1)=AMN
         LB(I2)=13
         E(I2)=DM
          ELSE
          LB(I2)=2
          E(I2)=AMN
         LB(I1)=13
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
          IF(M12.EQ.43)THEN
          IF(iabs(LB(I1)).EQ.8)THEN
          LB(I1)=1
          E(I1)=AMP
         LB(I2)=10
         E(I2)=DM
          ELSE
          LB(I2)=1
          E(I2)=AMP
         LB(I1)=10
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
          IF(M12.EQ.44)THEN
          IF(iabs(LB(I1)).EQ.8)THEN
          LB(I1)=1
          E(I1)=AMP
         LB(I2)=12
         E(I2)=DM
          ELSE
          LB(I2)=1
          E(I2)=AMP
         LB(I1)=12
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
          IF(M12.EQ.46)THEN
          IF(iabs(LB(I1)).EQ.6)THEN
          LB(I1)=2
          E(I1)=AMN
         LB(I2)=10
         E(I2)=DM
          ELSE
          LB(I2)=2
          E(I2)=AMN
         LB(I1)=10
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
          IF(M12.EQ.47)THEN
          IF(iabs(LB(I1)).EQ.6)THEN
          LB(I1)=2
          E(I1)=AMN
         LB(I2)=12
         E(I2)=DM
          ELSE
          LB(I2)=2
          E(I2)=AMN
         LB(I1)=12
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
          IF(M12.EQ.48)THEN
          IF(iabs(LB(I1)).EQ.7)THEN
          LB(I1)=2
          E(I1)=AMN
         LB(I2)=11
         E(I2)=DM
          ELSE
          LB(I2)=2
          E(I2)=AMN
         LB(I1)=11
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
          IF(M12.EQ.49)THEN
          IF(iabs(LB(I1)).EQ.7)THEN
          LB(I1)=2
          E(I1)=AMN
         LB(I2)=12
         E(I2)=DM
          ELSE
          LB(I2)=2
          E(I2)=AMN
         LB(I1)=12
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
          IF(M12.EQ.50)THEN
          IF(iabs(LB(I1)).EQ.7)THEN
          LB(I1)=1
          E(I1)=AMP
         LB(I2)=10
         E(I2)=DM
          ELSE
          LB(I2)=1
          E(I2)=AMP
         LB(I1)=10
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
          IF(M12.EQ.51)THEN
          IF(iabs(LB(I1)).EQ.7)THEN
          LB(I1)=2
          E(I1)=AMN
         LB(I2)=11
         E(I2)=DM
          ELSE
          LB(I2)=2
          E(I2)=AMN
         LB(I1)=11
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
          IF(M12.EQ.52)THEN
          IF(iabs(LB(I1)).EQ.7)THEN
          LB(I1)=1
          E(I1)=AMP
         LB(I2)=12
         E(I2)=DM
          ELSE
          LB(I2)=1
          E(I2)=AMP
         LB(I1)=12
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
          IF(M12.EQ.53)THEN
          IF(iabs(LB(I1)).EQ.7)THEN
          LB(I1)=2
          E(I1)=AMN
         LB(I2)=13
         E(I2)=DM
          ELSE
          LB(I2)=2
          E(I2)=AMN
         LB(I1)=13
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
          IF(M12.EQ.54)THEN
          IF(iabs(LB(I1)).EQ.10)THEN
          LB(I1)=2
          E(I1)=AMN
         LB(I2)=13
         E(I2)=DM
          ELSE
          LB(I2)=2
          E(I2)=AMN
         LB(I1)=13
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
          IF(M12.EQ.55)THEN
          IF(iabs(LB(I1)).EQ.10)THEN
          LB(I1)=1
          E(I1)=AMP
         LB(I2)=12
         E(I2)=DM
          ELSE
          LB(I2)=1
          E(I2)=AMP
         LB(I1)=12
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
          IF(M12.EQ.56)THEN
          IF(iabs(LB(I1)).EQ.11)THEN
          LB(I1)=2
          E(I1)=AMN
         LB(I2)=13
         E(I2)=DM
          ELSE
          LB(I2)=2
          E(I2)=AMN
         LB(I1)=13
         E(I1)=DM
          ENDIF
         GO TO 207
         ENDIF
          IF(M12.EQ.57)THEN
          IF(iabs(LB(I1)).EQ.11)THEN
          LB(I1)=1
          E(I1)=AMP
         LB(I2)=12
         E(I2)=DM
          ELSE
          LB(I2)=1
          E(I2)=AMP
         LB(I1)=12
         E(I1)=DM
          ENDIF
         ENDIF
          GO TO 207
206       IF(M12.EQ.1)THEN
          IF(iabs(LB(I1)).EQ.8)THEN
          LB(I2)=2
          LB(I1)=1
          E(I1)=AMP
          ELSE
          LB(I1)=2
          LB(I2)=1
          E(I2)=AMP
          ENDIF
         GO TO 207
          ENDIF
          IF(M12.EQ.2)THEN
          IF(iabs(LB(I1)).EQ.7)THEN
          LB(I2)=1
          LB(I1)=2
          E(I1)=AMN
          ELSE
          LB(I1)=1
          LB(I2)=2
          E(I2)=AMN
          ENDIF
         GO TO 207
          ENDIF
          IF(M12.EQ.3)THEN
          LB(I1)=1
          LB(I2)=1
          E(I1)=AMP
          E(I2)=AMP
         GO TO 207
          ENDIF
          IF(M12.EQ.4)THEN
          LB(I1)=1
          LB(I2)=1
          E(I1)=AMP
          E(I2)=AMP
         GO TO 207
          ENDIF
          IF(M12.EQ.5)THEN
          LB(I1)=2
          LB(I2)=2
          E(I1)=AMN
          E(I2)=AMN
         GO TO 207
          ENDIF
          IF(M12.EQ.6)THEN
          LB(I1)=2
          LB(I2)=2
          E(I1)=AMN
          E(I2)=AMN
         GO TO 207
          ENDIF
          IF(M12.EQ.7)THEN
          IF(iabs(LB(I1)).EQ.1)THEN
          LB(I1)=1
          LB(I2)=2
          E(I1)=AMP
          E(I2)=AMN
          ELSE
          LB(I1)=2
          LB(I2)=1
          E(I1)=AMN
          E(I2)=AMP
          ENDIF
         GO TO 207
          ENDIF
          IF(M12.EQ.8)THEN
          IF(iabs(LB(I1)).EQ.2)THEN
          LB(I1)=2
          LB(I2)=1
          E(I1)=AMN
          E(I2)=AMP
          ELSE
          LB(I1)=1
          LB(I2)=2
          E(I1)=AMP
          E(I2)=AMN
          ENDIF
         GO TO 207
          ENDIF
          IF(M12.EQ.9)THEN
          LB(I1)=1
          LB(I2)=1
          E(I1)=AMP
          E(I2)=AMP
         GO TO 207
         ENDIF
          IF(M12.EQ.12)THEN
          LB(I1)=2
          LB(I2)=1
          E(I1)=AMN
          E(I2)=AMP
         GO TO 207
         ENDIF
          IF(M12.EQ.11)THEN
          LB(I1)=2
          LB(I2)=1
          E(I1)=AMN
          E(I2)=AMP
         GO TO 207
         ENDIF
          IF(M12.EQ.12)THEN
          LB(I1)=1
          LB(I2)=2
          E(I1)=AMP
          E(I2)=AMN
         ENDIF
207       PR   = PRF
          C1   = 1.0 - 2.0 * RANART(NSEED)
              if(srt.le.2.14)C1= 1.0 - 2.0 * RANART(NSEED)
         if(srt.gt.2.14.and.srt.le.2.4)c1=ang(srt,iseed)
         if(srt.gt.2.4)then
             xptr=0.33*pr
         cc1=ptr(xptr,iseed)
         scheck=pr**2-cc1**2
         if(scheck.lt.0) then
            write(99,*) 'scheck4: ', scheck
            scheck=0.
         endif
         c1=sqrt(scheck)/pr
         endif
          T1   = 2.0 * PI * RANART(NSEED)
          IBLOCK=3
      ENDIF
      if(ianti.eq.1 .and. lb(i1).ge.1 .and. lb(i2).ge.1)then
         lb(i1) = -lb(i1)
         lb(i2) = -lb(i2)
      endif
 107  IF(PX .EQ. 0.0 .AND. PY .EQ. 0.0) THEN
         T2 = 0.0
      ELSE
         T2=ATAN2(PY,PX)
      END IF
      scheck=1.0 - C1**2
      if(scheck.lt.0) then
         write(99,*) 'scheck5: ', scheck
         scheck=0.
      endif
      S1=SQRT(scheck)
      scheck=1.0 - C2**2
      if(scheck.lt.0) then
         write(99,*) 'scheck6: ', scheck
         scheck=0.
      endif
      S2=SQRT(scheck)
      CT1  = COS(T1)
      ST1  = SIN(T1)
      CT2  = COS(T2)
      ST2  = SIN(T2)
      PZ   = PR * ( C1*C2 - S1*S2*CT1 )
      SS   = C2 * S1 * CT1  +  S2 * C1
      PX   = PR * ( SS*CT2 - S1*ST1*ST2 )
      PY   = PR * ( SS*ST2 + S1*ST1*CT2 )
      RETURN
306     CONTINUE
              if(XSK5/sigK.gt.RANART(NSEED))then
              pz1=p(3,i1)
              pz2=p(3,i2)
                LB(I1) = 1 + int(2 * RANART(NSEED))
                LB(I2) = 1 + int(2 * RANART(NSEED))
              nnn=nnn+1
                LPION(NNN,IRUN)=29
                EPION(NNN,IRUN)=APHI
                iblock = 222
              GO TO 208
               ENDIF
                IBLOCK=11
                if(ianti .eq. 1)iblock=-11
              pz1=p(3,i1)
              pz2=p(3,i2)
              nnn=nnn+1
                LPION(NNN,IRUN)=23
                EPION(NNN,IRUN)=Aka
              if(srt.le.2.63)then
              ic=1
                LB(I1) = 1 + int(2 * RANART(NSEED))
              LB(I2)=14
              GO TO 208
                ENDIF
       if(srt.le.2.74.and.srt.gt.2.63)then
              if(XSK1/(XSK1+XSK2).gt.RANART(NSEED))then
              ic=1
                LB(I1) = 1 + int(2 * RANART(NSEED))
              LB(I2)=14
              else
                   LB(I1) = 1 + int(2 * RANART(NSEED))
                   LB(I2) = 15 + int(3 * RANART(NSEED))
              ic=2
              endif
              GO TO 208
       endif
       if(srt.le.2.77.and.srt.gt.2.74)then
              if(xsk1/(xsk1+xsk2+xsk3).
     1          gt.RANART(NSEED))then
              ic=1
                LB(I1) = 1 + int(2 * RANART(NSEED))
              LB(I2)=14
              go to 208
              else
              if(xsk2/(xsk2+xsk3).gt.RANART(NSEED))then
              ic=2
                LB(I1) = 1 + int(2 * RANART(NSEED))
                LB(I2) = 15 + int(3 * RANART(NSEED))
              else
              ic=3
                LB(I1) = 6 + int(4 * RANART(NSEED))
              lb(i2)=14
              endif
              GO TO 208
              endif
       endif
       if(srt.gt.2.77)then
              if(xsk1/(xsk1+xsk2+xsk3+xsk4).gt.RANART(NSEED))then
              ic=1
                LB(I1) = 1 + int(2 * RANART(NSEED))
              LB(I2)=14
              go to 208
       else
          if(xsk3/(xsk2+xsk3+xsk4).gt.RANART(NSEED))then
              ic=3
                LB(I1) = 6 + int(4 * RANART(NSEED))
              lb(i2)=14
              go to 208
          else
              if(xsk2/(xsk2+xsk4).gt.RANART(NSEED))then
                   LB(I1) = 1 + int(2 * RANART(NSEED))
                   LB(I2) = 15 + int(3 * RANART(NSEED))
              ic=2
              else
              ic=4
                LB(I1) = 6 + int(4 * RANART(NSEED))
                LB(I2) = 15 + int(3 * RANART(NSEED))
              endif
              go to 208
          endif
       endif
       endif
208             continue
         if(ianti.eq.1 .and. lb(i1).ge.1 .and. lb(i2).ge.1)then
          lb(i1) = - lb(i1)
          lb(i2) = - lb(i2)
          if(LPION(NNN,IRUN) .eq. 23)LPION(NNN,IRUN)=21
         endif
       lbi1=lb(i1)
       lbi2=lb(i2)
           NTRY1=0
128        CALL BBKAON(ic,SRT,PX3,PY3,PZ3,DM3,PX4,PY4,PZ4,DM4,
     &  PPX,PPY,PPZ,icou1)
       NTRY1=NTRY1+1
       if((icou1.lt.0).AND.(NTRY1.LE.20))GO TO 128
       CALL ROTATE(PX,PY,PZ,PX3,PY3,PZ3)
       CALL ROTATE(PX,PY,PZ,PX4,PY4,PZ4)
       CALL ROTATE(PX,PY,PZ,PPX,PPY,PPZ)
              E1CM    = SQRT (dm3**2 + PX3**2 + PY3**2 + PZ3**2)
              P1BETA  = PX3*BETAX + PY3*BETAY + PZ3*BETAZ
              TRANSF  = GAMMA * ( GAMMA * P1BETA / (GAMMA + 1) + E1CM )
              Pt1i1 = BETAX * TRANSF + PX3
              Pt2i1 = BETAY * TRANSF + PY3
              Pt3i1 = BETAZ * TRANSF + PZ3
             Eti1   = DM3
                E2CM    = SQRT (dm4**2 + PX4**2 + PY4**2 + PZ4**2)
                P2BETA  = PX4*BETAX+PY4*BETAY+PZ4*BETAZ
                TRANSF  = GAMMA * (GAMMA*P2BETA / (GAMMA + 1.) + E2CM)
                Pt1I2 = BETAX * TRANSF + PX4
                Pt2I2 = BETAY * TRANSF + PY4
                Pt3I2 = BETAZ * TRANSF + PZ4
              EtI2   = DM4
                EPCM=SQRT(aka**2+PPX**2+PPY**2+PPZ**2)
                PPBETA=PPX*BETAX+PPY*BETAY+PPZ*BETAZ
                TRANSF=GAMMA*(GAMMA*PPBETA/(GAMMA+1.)+EPCM)
                PPION(1,NNN,IRUN)=BETAX*TRANSF+PPX
                PPION(2,NNN,IRUN)=BETAY*TRANSF+PPY
                PPION(3,NNN,IRUN)=BETAZ*TRANSF+PPZ
                dppion(nnn,irun)=dpertp(i1)*dpertp(i2)
                    RPION(1,NNN,IRUN)=R(1,I1)
                    RPION(2,NNN,IRUN)=R(2,I1)
                    RPION(3,NNN,IRUN)=R(3,I1)
              p(1,i1)=pt1i1
              p(2,i1)=pt2i1
              p(3,i1)=pt3i1
              e(i1)=eti1
              lb(i1)=lbi1
              p(1,i2)=pt1i2
              p(2,i2)=pt2i2
              p(3,i2)=pt3i2
              e(i2)=eti2
              lb(i2)=lbi2
                PX1     = P(1,I1)
                PY1     = P(2,I1)
                PZ1     = P(3,I1)
              EM1       = E(I1)
                ID(I1)  = 2
                ID(I2)  = 2
                ID1     = ID(I1)
                if(LPION(NNN,IRUN) .ne. 29) IBLOCK=11
        LB1=LB(I1)
        LB2=LB(I2)
        AM1=EM1
       am2=em2
        E1= SQRT( EM1**2 + PX1**2 + PY1**2 + PZ1**2 )
       RETURN
 108   CONTINUE
           if(idpert.eq.1.and.ipert1.eq.1.and.npertd.ge.1) then
              ndloop=npertd
           elseif(idpert.eq.2.and.npertd.ge.1) then
              ndloop=npertd+1
           else
              ndloop=1
           endif
           dprob1=sdprod/sig/float(npertd)
           do idloop=1,ndloop
              CALL bbdangle(pxd,pyd,pzd,nt,ipert1,ianti,idloop,pfinal,
     1 dprob1,lbm)
              CALL ROTATE(PX,PY,PZ,PXd,PYd,PZd)
              xmass=xmd
              E1dCM=SQRT(xmass**2+PXd**2+PYd**2+PZd**2)
              P1dBETA=PXd*BETAX+PYd*BETAY+PZd*BETAZ
              TRANSF=GAMMA*(GAMMA*P1dBETA/(GAMMA+1.)+E1dCM)
              pxi1=BETAX*TRANSF+PXd
              pyi1=BETAY*TRANSF+PYd
              pzi1=BETAZ*TRANSF+PZd
              if(ianti.eq.0)then
                 lbd=42
              else
                 lbd=-42
              endif
              if(idpert.eq.1.and.ipert1.eq.1.and.npertd.ge.1) then
                 nnn=nnn+1
                 PPION(1,NNN,IRUN)=pxi1
                 PPION(2,NNN,IRUN)=pyi1
                 PPION(3,NNN,IRUN)=pzi1
                 EPION(NNN,IRUN)=xmd
                 LPION(NNN,IRUN)=lbd
                 RPION(1,NNN,IRUN)=R(1,I1)
                 RPION(2,NNN,IRUN)=R(2,I1)
                 RPION(3,NNN,IRUN)=R(3,I1)
                 dppion(NNN,IRUN)=sdprod/sig/float(npertd)
              elseif(idpert.eq.2.and.idloop.le.npertd) then
                 ppd(1,idloop)=pxi1
                 ppd(2,idloop)=pyi1
                 ppd(3,idloop)=pzi1
                 lbpd(idloop)=lbd
              else
                 E(i1)=xmm
                 E2piCM=SQRT(xmm**2+PXd**2+PYd**2+PZd**2)
                 P2piBETA=-PXd*BETAX-PYd*BETAY-PZd*BETAZ
                 TRANSF=GAMMA*(GAMMA*P2piBETA/(GAMMA+1.)+E2piCM)
                 pxi2=BETAX*TRANSF-PXd
                 pyi2=BETAY*TRANSF-PYd
                 pzi2=BETAZ*TRANSF-PZd
                 p(1,i1)=pxi2
                 p(2,i1)=pyi2
                 p(3,i1)=pzi2
                 LB(I1)=lbm
                 PX1=P(1,I1)
                 PY1=P(2,I1)
                 PZ1=P(3,I1)
                 EM1=E(I1)
                 ID(I1)=2
                 ID1=ID(I1)
                 E1=SQRT(EM1**2+PX1**2+PY1**2+PZ1**2)
                 lb1=lb(i1)
                 p(1,i2)=pxi1
                 p(2,i2)=pyi1
                 p(3,i2)=pzi1
                 lb(i2)=lbd
                 lb2=lb(i2)
                 E(i2)=xmd
                 EtI2=E(I2)
                 ID(I2)=2
                 if(idpert.eq.2.and.idloop.eq.ndloop) then
                    do ipertd=1,npertd
                       nnn=nnn+1
                       PPION(1,NNN,IRUN)=ppd(1,ipertd)
                       PPION(2,NNN,IRUN)=ppd(2,ipertd)
                       PPION(3,NNN,IRUN)=ppd(3,ipertd)
                       EPION(NNN,IRUN)=xmd
                       LPION(NNN,IRUN)=lbpd(ipertd)
                       RPION(1,NNN,IRUN)=R(1,I1)
                       RPION(2,NNN,IRUN)=R(2,I1)
                       RPION(3,NNN,IRUN)=R(3,I1)
                       dppion(NNN,IRUN)=1./float(npertd)
                    enddo
                 endif
              endif
           enddo
           IBLOCK=501
           return
      END
