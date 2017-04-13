      SUBROUTINE CRNN(IRUN,PX,PY,PZ,SRT,I1,I2,IBLOCK,
     1NTAG,SIGNN,SIG,NT,ipert1)
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,aka=0.498,AP2=0.13957,AM0=1.232,
     2  PI=3.1415926,CUTOFF=1.8966,AVMASS=0.9383,APHI=1.020)
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
        COMMON/TABLE/ xarray(0:1000),earray(0:1000)
        common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      common/leadng/lb1,px1,py1,pz1,em1,e1,xfnl,yfnl,zfnl,tfnl,
     1 px1n,py1n,pz1n,dp1n
      COMMON/RNDF77/NSEED
      common /dpi/em2,lb2
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
      common /para8/ idpert,npertd,idxsec
      dimension ppd(3,npdmax),lbpd(npdmax)
      SAVE   
      n12=0
      m12=0
      IBLOCK=0
      NTAG=0
      EM1=E(I1)
      EM2=E(I2)
      PR=SQRT( PX**2 + PY**2 + PZ**2 )
      C2=PZ / PR
      X1=RANART(NSEED)
      ianti=0
      if(lb(i1).lt.0 .and. lb(i2).lt.0) ianti=1
      call sbbdm(srt,sdprod,ianti,lbm,xmm,pfinal)
      if(idpert.eq.1.and.ipert1.eq.1) then
         IF (SRT .LT. 2.012) RETURN
         if((iabs(lb(i1)).eq.1.or.iabs(lb(i1)).eq.2)
     1        .and.(iabs(lb(i2)).eq.1.or.iabs(lb(i2)).eq.2)) then
            goto 108
         else
            return
         endif
      endif
      IF (X1.LE.(SIGNN/SIG)) THEN
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
         IF (SRT .LT. 2.012) RETURN
       call N1535(iabs(lb(i1)),iabs(lb(i2)),srt,x1535)
       SIG3=3.*(X3pi(SRT)+x33pi(srt))
       SIG4=4.*X2pi(srt)
       s4pi=x4pi(srt)
       srho=xrho(srt)
       somega=omega(srt)
       akp=0.498
       ak0=0.498
       ana=0.94
       ada=1.232
       al=1.1157
       as=1.1197
       xsk1=0
       xsk2=0
       xsk3=0
       xsk4=0
       xsk5=0
       t1nlk=ana+al+akp
       if(srt.le.t1nlk)go to 222
       XSK1=1.5*PPLPK(SRT)
       t1dlk=ada+al+akp
       t2dlk=ada+al-akp
       if(srt.le.t1dlk)go to 222
       es=srt
       pmdlk2=(es**2-t1dlk**2)*(es**2-t2dlk**2)/(4.*es**2)
       pmdlk=sqrt(pmdlk2)
       XSK3=1.5*PPLPK(srt)
       t1nsk=ana+as+akp
       t2nsk=ana+as-akp
       if(srt.le.t1nsk)go to 222
       pmnsk2=(es**2-t1nsk**2)*(es**2-t2nsk**2)/(4.*es**2)
       pmnsk=sqrt(pmnsk2)
       XSK2=1.5*(PPK1(srt)+PPK0(srt))
       t1DSk=aDa+aS+akp
       t2DSk=aDa+aS-akp
       if(srt.le.t1dsk)go to 222
       pmDSk2=(es**2-t1DSk**2)*(es**2-t2DSk**2)/(4.*es**2)
       pmDSk=sqrt(pmDSk2)
       XSK4=1.5*(PPK1(srt)+PPK0(srt))
       if(srt.le.(2.*amn+aphi))go to 222
       xsk5 = 0.0001
 222   SIGK=XSK1+XSK2+XSK3+XSK4
        XSK1 = 2.0 * XSK1
        XSK2 = 2.0 * XSK2
        XSK3 = 2.0 * XSK3
        XSK4 = 2.0 * XSK4
        SIGK = 2.0 * SIGK + xsk5
        lb1=iabs(lb(i1))
        lb2=iabs(lb(i2))
        IF((LB(I1)*LB(I2).EQ.1).or.
     &       ((lb1.le.17.and.lb1.ge.14).and.(lb2.le.17.and.lb2.ge.14)).
     &       or.((lb1.le.2).and.(lb2.le.17.and.lb2.ge.14)).
     &       or.((lb2.le.2).and.(lb1.le.17.and.lb1.ge.14)))THEN
           IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
           SIG1=SIGMA(SRT,1,1,0)+0.5*SIGMA(SRT,1,1,1)
           SIG2=1.5*SIGMA(SRT,1,1,1)
           SIGND=SIG1+SIG2+SIG3+SIG4+X1535+SIGK+s4pi+srho+somega
           IF (X1.GT.(SIGNN+SIGND+sdprod)/SIG)RETURN
           DIR=SIG3/SIGND
           IF(RANART(NSEED).LE.DIR)GO TO 106
           IF(RANART(NSEED).LE.SIGK/(SIGK+X1535+SIG4+SIG2+SIG1
     &          +s4pi+srho+somega))GO TO 306
           if(RANART(NSEED).le.s4pi/(x1535+sig4+sig2+sig1
     &          +s4pi+srho+somega))go to 307
           if(RANART(NSEED).le.srho/(x1535+sig4+sig2+sig1
     &          +srho+somega))go to 308
           if(RANART(NSEED).le.somega/(x1535+sig4+sig2+sig1
     &          +somega))go to 309
           if(RANART(NSEED).le.x1535/(sig1+sig2+sig4+x1535))then
              N12=9
           ELSE 
              IF(RANART(NSEED).LE.SIG4/(SIG1+sig2+sig4))THEN
                 N12=66
                 GO TO 1012
              else
                 N12=3
                 IF (RANART(NSEED).GT.SIG1/(SIG1+SIG2))N12=4
              ENDIF
           endif
           GO TO 1011
        ENDIF
        IF(iabs(LB(I1)).EQ.2.AND.iabs(LB(I2)).EQ.2)THEN
           IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
           SIG1=SIGMA(SRT,1,1,0)+0.5*SIGMA(SRT,1,1,1)
           SIG2=1.5*SIGMA(SRT,1,1,1)
           SIGND=SIG1+SIG2+X1535+SIG3+SIG4+SIGK+s4pi+srho+somega
           IF (X1.GT.(SIGNN+SIGND+sdprod)/SIG)RETURN
           dir=sig3/signd
           IF(RANART(NSEED).LE.DIR)GO TO 106
           IF(RANART(NSEED).LE.SIGK/(SIGK+X1535+SIG4+SIG2+SIG1
     &          +s4pi+srho+somega))GO TO 306
           if(RANART(NSEED).le.s4pi/(x1535+sig4+sig2+sig1
     &          +s4pi+srho+somega))go to 307
           if(RANART(NSEED).le.srho/(x1535+sig4+sig2+sig1
     &          +srho+somega))go to 308
           if(RANART(NSEED).le.somega/(x1535+sig4+sig2+sig1
     &          +somega))go to 309
           IF(RANART(NSEED).LE.X1535/(x1535+sig1+sig2+sig4))THEN
              N12=10
           ELSE 
              if(RANART(NSEED).le.sig4/(sig1+sig2+sig4))then
                 N12=67
                 GO TO 1013
              else
                 N12=6
                 IF (RANART(NSEED).GT.SIG1/(SIG1+SIG2))N12=5
              ENDIF
           endif
           GO TO 1011
        ENDIF
        IF(LB(I1)*LB(I2).EQ.2)THEN
           IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
           SIG1=0.5*SIGMA(SRT,1,1,1)+0.25*SIGMA(SRT,1,1,0)
           IF(NSTAR.EQ.1)THEN
              SIG2=(3./4.)*SIGMA(SRT,2,0,1)
           ELSE
              SIG2=0.
           ENDIF
           SIGND=2.*(SIG1+SIG2+X1535)+sig3+sig4+SIGK+s4pi+srho+somega
           IF (X1.GT.(SIGNN+SIGND+sdprod)/SIG)RETURN
           dir=sig3/signd
           IF(RANART(NSEED).LE.DIR)GO TO 106
           IF(RANART(NSEED).LE.SIGK/(SIGND-SIG3))GO TO 306
           if(RANART(NSEED).le.s4pi/(signd-sig3-sigk))go to 307
           if(RANART(NSEED).le.srho/(signd-sig3-sigk-s4pi))go to 308
           if(RANART(NSEED).le.somega/(signd-sig3-sigk-s4pi-srho))
     1          go to 309
           IF(RANART(NSEED).LT.X1535/(SIG1+SIG2+X1535+0.5*sig4))THEN
              N12=11
              IF(RANART(NSEED).LE.0.5)N12=12
           ELSE 
              if(RANART(NSEED).le.sig4/(sig4+2.*(sig1+sig2)))then
                 N12=68
                 GO TO 1014
              else
                 IF(RANART(NSEED).LE.SIG1/(SIG1+SIG2))THEN
                    N12=2
                    IF(RANART(NSEED).GE.0.5)N12=1
                 ELSE
                    N12=8
                    IF(RANART(NSEED).GE.0.5)N12=7
                 ENDIF
              ENDIF
           ENDIF
        endif
 1011   iblock=2
        CONTINUE
          DMAX = SRT - AVMASS-0.005
          DMAX = SRT - AVMASS-0.005
          DMIN = 1.078
                   IF(N12.LT.7)THEN
          IF(DMAX.LT.1.232) THEN
          FM=FDE(DMAX,SRT,0.)
          ELSE
             xdmass=1.232
          FM=FDE(xdmass,SRT,1.)
          ENDIF
          IF(FM.EQ.0.)FM=1.E-09
          NTRY1=0
10        DM = RANART(NSEED) * (DMAX-DMIN) + DMIN
          NTRY1=NTRY1+1
          IF((RANART(NSEED) .GT. FDE(DM,SRT,1.)/FM).AND.
     1    (NTRY1.LE.30)) GOTO 10
          if(dm.gt.1.47) goto 10
              GO TO 13
              ENDIF
                   IF((n12.eq.7).or.(n12.eq.8))THEN
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
              ENDIF
                    IF(n12.ge.17)then
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
         GO TO 13
             ENDIF
1012       iblock=43
       call Rmasdd(srt,1.232,1.232,1.08,
     &  1.08,ISEED,1,dm1,dm2)
       call Rmasdd(srt,1.232,1.44,1.08,
     &  1.08,ISEED,3,dm1n,dm2n)
       IF(N12.EQ.66)THEN
       XFINAL=RANART(NSEED)
       IF(XFINAL.LE.0.25)THEN
       LB(I1)=9
       LB(I2)=7
       e(i1)=dm1
       e(i2)=dm2
       GO TO 200
       ENDIF
       IF((XFINAL.gt.0.25).and.(xfinal.le.0.5))THEN
       LB(I1)=8
       LB(I2)=8
       e(i1)=dm1
       e(i2)=dm2
       GO TO 200
       ENDIF
       IF((XFINAL.gt.0.5).and.(xfinal.le.0.75))THEN
       LB(I1)=9
       LB(I2)=10
       e(i1)=dm1n
       e(i2)=dm2n
       GO TO 200
       ENDIF
       IF(XFINAL.gt.0.75)then
       LB(I1)=8
       LB(I2)=11
       e(i1)=dm1n
       e(i2)=dm2n
       GO TO 200
       ENDIF
       ENDIF
1013       iblock=43
       call Rmasdd(srt,1.232,1.232,1.08,
     &  1.08,ISEED,1,dm1,dm2)
       call Rmasdd(srt,1.232,1.44,1.08,
     &  1.08,ISEED,3,dm1n,dm2n)
       IF(N12.EQ.67)THEN
       XFINAL=RANART(NSEED)
       IF(XFINAL.LE.0.25)THEN
       LB(I1)=7
       LB(I2)=7
       e(i1)=dm1
       e(i2)=dm2
       GO TO 200
        ENDIF
       IF((XFINAL.gt.0.25).and.(xfinal.le.0.5))THEN
       LB(I1)=6
       LB(I2)=8
       e(i1)=dm1
       e(i2)=dm2
       GO TO 200
       ENDIF
       IF((XFINAL.gt.0.5).and.(xfinal.le.0.75))THEN
       LB(I1)=7
       LB(I2)=10
       e(i1)=dm1n
       e(i2)=dm2n
       GO TO 200
       ENDIF
       IF(XFINAL.gt.0.75)then
       LB(I1)=8
       LB(I2)=11
       e(i1)=dm1n
       e(i2)=dm2n
       GO TO 200
       ENDIF
       ENDIF
1014       iblock=43
       call Rmasdd(srt,1.232,1.232,1.08,
     &  1.08,ISEED,1,dm1,dm2)
       call Rmasdd(srt,1.232,1.44,1.08,
     &  1.08,ISEED,3,dm1n,dm2n)
       IF(N12.EQ.68)THEN
       XFINAL=RANART(NSEED)
       IF(XFINAL.LE.0.25)THEN
       LB(I1)=7
       LB(I2)=8
       e(i1)=dm1
       e(i2)=dm2
       GO TO 200
       ENDIF
       IF((XFINAL.gt.0.25).and.(xfinal.le.0.5))THEN
       LB(I1)=9
       LB(I2)=6
       e(i1)=dm1
       e(i2)=dm2
       GO TO 200
       ENDIF
       IF((XFINAL.gt.0.5).and.(xfinal.le.0.75))THEN
       LB(I1)=7
       LB(I2)=11
       e(i1)=dm1n
       e(i2)=dm2n
       GO TO 200
       ENDIF
       IF(XFINAL.gt.0.75)then
       LB(I1)=8
       LB(I2)=10
       e(i1)=dm1n
       e(i2)=dm2n
       GO TO 200
       ENDIF
       ENDIF
13       CONTINUE
          IF(N12.EQ.1)THEN
          IF(iabs(LB(I1)).EQ.1)THEN
          LB(I2)=2
          LB(I1)=8
          E(I1)=DM
          ELSE
          LB(I1)=2
          LB(I2)=8
          E(I2)=DM
          ENDIF
         GO TO 200
          ENDIF
          IF(N12.EQ.2)THEN
          IF(iabs(LB(I1)).EQ.2)THEN
          LB(I2)=1
          LB(I1)=7
          E(I1)=DM
          ELSE
          LB(I1)=1
          LB(I2)=7
          E(I2)=DM
          ENDIF
         GO TO 200
          ENDIF
          IF(N12.EQ.3)THEN
          LB(I1)=9
          E(I1)=DM
          LB(I2)=2
          E(I2)=AMN
         GO TO 200
          ENDIF
          IF(N12.EQ.4)THEN
          LB(I2)=1
          LB(I1)=8
          E(I1)=DM
         GO TO 200
          ENDIF
          IF(N12.EQ.5)THEN
          LB(I2)=2
          LB(I1)=7
          E(I1)=DM
         GO TO 200
          ENDIF
          IF(N12.EQ.6)THEN
          LB(I1)=6
          E(I1)=DM
          LB(I2)=1
          E(I2)=AMP
         GO TO 200
          ENDIF
          IF(N12.EQ.7)THEN
          IF(iabs(LB(I1)).EQ.1)THEN
          LB(I1)=1
          LB(I2)=10
          E(I2)=DM
          ELSE
          LB(I2)=1
          LB(I1)=10
          E(I1)=DM
          ENDIF
         GO TO 200
          ENDIF
          IF(N12.EQ.8)THEN
          IF(iabs(LB(I1)).EQ.1)THEN
          LB(I2)=2
          LB(I1)=11
          E(I1)=DM
          ELSE
          LB(I1)=2
          LB(I2)=11
          E(I2)=DM
          ENDIF
         GO TO 200
          ENDIF
          IF(N12.EQ.9)THEN
          IF(RANART(NSEED).le.0.5)THEN
          LB(I2)=1
          LB(I1)=13
          E(I1)=DM
          ELSE
          LB(I1)=1
          LB(I2)=13
          E(I2)=DM
          ENDIF
         GO TO 200
          ENDIF
          IF(N12.EQ.10)THEN
          IF(RANART(NSEED).le.0.5)THEN
          LB(I2)=2
          LB(I1)=12
          E(I1)=DM
          ELSE
          LB(I1)=2
          LB(I2)=12
          E(I2)=DM
          ENDIF
         GO TO 200
          ENDIF
          IF(N12.EQ.11)THEN
          IF(iabs(LB(I1)).EQ.2)THEN
          LB(I1)=2
          LB(I2)=13
          E(I2)=DM
          ELSE
          LB(I2)=2
          LB(I1)=13
          E(I1)=DM
          ENDIF
         GO TO 200
          ENDIF
          IF(N12.EQ.12)THEN
          IF(iabs(LB(I1)).EQ.1)THEN
          LB(I1)=1
          LB(I2)=12
          E(I2)=DM
          ELSE
          LB(I2)=1
          LB(I1)=12
          E(I1)=DM
          ENDIF
          ENDIF
         endif
200       EM1=E(I1)
          EM2=E(I2)
          PR2   = (SRT**2 - EM1**2 - EM2**2)**2
     1                - 4.0 * (EM1*EM2)**2
          IF(PR2.LE.0.)PR2=1.e-09
          PR=SQRT(PR2)/(2.*SRT)
              if(srt.le.2.14)C1= 1.0 - 2.0 * RANART(NSEED)
         if(srt.gt.2.14.and.srt.le.2.4)c1=ang(srt,iseed)
         if(srt.gt.2.4)then
             xptr=0.33*pr
             cc1=ptr(xptr,iseed)
             scheck=pr**2-cc1**2
             if(scheck.lt.0) then
                write(99,*) 'scheck2: ', scheck
                scheck=0.
             endif
             c1=sqrt(scheck)/pr
         endif
          T1   = 2.0 * PI * RANART(NSEED)
       if(ianti.eq.1 .and. lb(i1).ge.1 .and. lb(i2).ge.1)then
         lb(i1) = -lb(i1)
         lb(i2) = -lb(i2)
       endif
          GO TO 107
106     CONTINUE
           NTRY1=0
123        CALL DDP2(SRT,ISEED,PX3,PY3,PZ3,DM3,PX4,PY4,PZ4,DM4,
     &  PPX,PPY,PPZ,icou1)
       NTRY1=NTRY1+1
       if((icou1.lt.0).AND.(NTRY1.LE.40))GO TO 123
       CALL ROTATE(PX,PY,PZ,PX3,PY3,PZ3)
       CALL ROTATE(PX,PY,PZ,PX4,PY4,PZ4)
       CALL ROTATE(PX,PY,PZ,PPX,PPY,PPZ)
                NNN=NNN+1
              XDIR=RANART(NSEED)
                IF(LB(I1)*LB(I2).EQ.1)THEN
                IF(XDIR.Le.0.2)then
                LPION(NNN,IRUN)=4
                EPION(NNN,IRUN)=AP1
              LB(I1)=9
              LB(I2)=7
       GO TO 205
                ENDIF
                IF((XDIR.LE.0.4).AND.(XDIR.GT.0.2))THEN
                LPION(NNN,IRUN)=4
                EPION(NNN,IRUN)=AP1
                LB(I1)=8
                LB(I2)=8
       GO TO 205
              ENDIF 
                IF((XDIR.LE.0.6).AND.(XDIR.GT.0.4))THEN
                LPION(NNN,IRUN)=3
                EPION(NNN,IRUN)=AP2
                LB(I1)=9
                LB(I2)=8
       GO TO 205
              ENDIF 
                IF((XDIR.LE.0.8).AND.(XDIR.GT.0.6))THEN
                LPION(NNN,IRUN)=5
                EPION(NNN,IRUN)=AP2
                LB(I1)=9
                LB(I2)=6
       GO TO 205
              ENDIF 
                IF(XDIR.GT.0.8)THEN
                LPION(NNN,IRUN)=5
                EPION(NNN,IRUN)=AP2
                LB(I1)=7
                LB(I2)=8
       GO TO 205
              ENDIF 
               ENDIF
                IF(iabs(LB(I1)).EQ.2.AND.iabs(LB(I2)).EQ.2)THEN
                IF(XDIR.Le.0.2)then
                LPION(NNN,IRUN)=4
                EPION(NNN,IRUN)=AP1
              LB(I1)=6
              LB(I2)=7
       GO TO 205
                ENDIF
                IF((XDIR.LE.0.4).AND.(XDIR.GT.0.2))THEN
                LPION(NNN,IRUN)=3
                EPION(NNN,IRUN)=AP2
                LB(I1)=6
                LB(I2)=9
       GO TO 205
              ENDIF 
                IF((XDIR.GT.0.4).AND.(XDIR.LE.0.6))THEN
                LPION(NNN,IRUN)=5
                EPION(NNN,IRUN)=AP2
                LB(I1)=9
                LB(I2)=8
       GO TO 205
              ENDIF 
                IF((XDIR.GT.0.6).AND.(XDIR.LE.0.8))THEN
                LPION(NNN,IRUN)=4
                EPION(NNN,IRUN)=AP1
                LB(I1)=7
                LB(I2)=7
       GO TO 205
              ENDIF 
                IF(XDIR.GT.0.8)THEN
                LPION(NNN,IRUN)=3
                EPION(NNN,IRUN)=AP2
                LB(I1)=7
                LB(I2)=8
       GO TO 205
              ENDIF 
              ENDIF
                IF(LB(I1)*LB(I2).EQ.2)THEN
                IF(XDIR.Le.0.17)then
                LPION(NNN,IRUN)=4
                EPION(NNN,IRUN)=AP1
              LB(I1)=6
              LB(I2)=9
       GO TO 205
                ENDIF
                IF((XDIR.LE.0.34).AND.(XDIR.GT.0.17))THEN
                LPION(NNN,IRUN)=3
                EPION(NNN,IRUN)=AP2
                LB(I1)=7
                LB(I2)=9
       GO TO 205
              ENDIF 
                IF((XDIR.GT.0.34).AND.(XDIR.LE.0.51))THEN
                LPION(NNN,IRUN)=5
                EPION(NNN,IRUN)=AP2
                LB(I1)=7
                LB(I2)=8
       GO TO 205
              ENDIF 
                IF((XDIR.GT.0.51).AND.(XDIR.LE.0.68))THEN
                LPION(NNN,IRUN)=3
                EPION(NNN,IRUN)=AP2
                LB(I1)=8
                LB(I2)=8
       GO TO 205
              ENDIF 
                IF((XDIR.GT.0.68).AND.(XDIR.LE.0.85))THEN
                LPION(NNN,IRUN)=4
                EPION(NNN,IRUN)=AP2
                LB(I1)=7
                LB(I2)=8
       GO TO 205
              ENDIF 
                IF(XDIR.GT.0.85)THEN
                LPION(NNN,IRUN)=5
                EPION(NNN,IRUN)=AP2
                LB(I1)=7
                LB(I2)=7
              ENDIF 
                ENDIF
205           E1CM    = SQRT (dm3**2 + PX3**2 + PY3**2 + PZ3**2)
              P1BETA  = PX3*BETAX + PY3*BETAY + PZ3*BETAZ
              TRANSF  = GAMMA * ( GAMMA * P1BETA / (GAMMA + 1) + E1CM )
              Pt1i1 = BETAX * TRANSF + PX3
              Pt2i1 = BETAY * TRANSF + PY3
              Pt3i1 = BETAZ * TRANSF + PZ3
             Eti1   = DM3
             if(ianti.eq.1 .and. lb(i1).ge.1 .and. lb(i2).ge.1)then
               lb(i1) = -lb(i1)
               lb(i2) = -lb(i2)
                if(LPION(NNN,IRUN) .eq. 3)then
                  LPION(NNN,IRUN)=5
                elseif(LPION(NNN,IRUN) .eq. 5)then
                  LPION(NNN,IRUN)=3
                endif
               endif
             lb1=lb(i1)
                E2CM    = SQRT (dm4**2 + PX4**2 + PY4**2 + PZ4**2)
                P2BETA  = PX4*BETAX+PY4*BETAY+PZ4*BETAZ
                TRANSF  = GAMMA * (GAMMA*P2BETA / (GAMMA + 1.) + E2CM)
                Pt1I2 = BETAX * TRANSF + PX4
                Pt2I2 = BETAY * TRANSF + PY4
                Pt3I2 = BETAZ * TRANSF + PZ4
              EtI2   = DM4
              lb2=lb(i2)
              p(1,i1)=pt1i1
              p(2,i1)=pt2i1
              p(3,i1)=pt3i1
              e(i1)=eti1
              lb(i1)=lb1
              p(1,i2)=pt1i2
              p(2,i2)=pt2i2
              p(3,i2)=pt3i2
              e(i2)=eti2
              lb(i2)=lb2
                PX1     = P(1,I1)
                PY1     = P(2,I1)
                PZ1     = P(3,I1)
              EM1       = E(I1)
                ID(I1)  = 2
                ID(I2)  = 2
                ID1     = ID(I1)
                IBLOCK=4
                EPCM=SQRT(EPION(NNN,IRUN)**2+PPX**2+PPY**2+PPZ**2)
                PPBETA=PPX*BETAX+PPY*BETAY+PPZ*BETAZ
                TRANSF=GAMMA*(GAMMA*PPBETA/(GAMMA+1.)+EPCM)
                PPION(1,NNN,IRUN)=BETAX*TRANSF+PPX
                PPION(2,NNN,IRUN)=BETAY*TRANSF+PPY
                PPION(3,NNN,IRUN)=BETAZ*TRANSF+PPZ
                dppion(nnn,irun)=dpertp(i1)*dpertp(i2)
                RPION(1,NNN,IRUN)=R(1,I1)
                RPION(2,NNN,IRUN)=R(2,I1)
                RPION(3,NNN,IRUN)=R(3,I1)
              go to 90005
 108       CONTINUE
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
           go to 90005
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
                 IBLOCK=9
                 if(ianti .eq. 1)iblock=-9
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
           NTRY1=0
127        CALL BBKAON(ic,SRT,PX3,PY3,PZ3,DM3,PX4,PY4,PZ4,DM4,
     &  PPX,PPY,PPZ,icou1)
       NTRY1=NTRY1+1
       if((icou1.lt.0).AND.(NTRY1.LE.20))GO TO 127
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
             lbi1=lb(i1)
                E2CM    = SQRT (dm4**2 + PX4**2 + PY4**2 + PZ4**2)
                P2BETA  = PX4*BETAX+PY4*BETAY+PZ4*BETAZ
                TRANSF  = GAMMA * (GAMMA*P2BETA / (GAMMA + 1.) + E2CM)
                Pt1I2 = BETAX * TRANSF + PX4
                Pt2I2 = BETAY * TRANSF + PY4
                Pt3I2 = BETAZ * TRANSF + PZ4
              EtI2   = DM4
              lbi2=lb(i2)
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
              go to 90005
307     CONTINUE
           NTRY1=0
125        CALL DDrho(SRT,ISEED,PX3,PY3,PZ3,DM3,PX4,PY4,PZ4,DM4,
     &  PPX,PPY,PPZ,amrho,icou1)
       NTRY1=NTRY1+1
       if((icou1.lt.0).AND.(NTRY1.LE.20))GO TO 125
       CALL ROTATE(PX,PY,PZ,PX3,PY3,PZ3)
       CALL ROTATE(PX,PY,PZ,PX4,PY4,PZ4)
       CALL ROTATE(PX,PY,PZ,PPX,PPY,PPZ)
                NNN=NNN+1
              arho=amrho
              XDIR=RANART(NSEED)
                IF(LB(I1)*LB(I2).EQ.1)THEN
                IF(XDIR.Le.0.2)then
                LPION(NNN,IRUN)=26
                EPION(NNN,IRUN)=Arho
              LB(I1)=9
              LB(I2)=7
       GO TO 2051
                ENDIF
                IF((XDIR.LE.0.4).AND.(XDIR.GT.0.2))THEN
                LPION(NNN,IRUN)=26
                EPION(NNN,IRUN)=Arho
                LB(I1)=8
                LB(I2)=8
       GO TO 2051
              ENDIF 
                IF((XDIR.LE.0.6).AND.(XDIR.GT.0.4))THEN
                LPION(NNN,IRUN)=25
                EPION(NNN,IRUN)=Arho
                LB(I1)=9
                LB(I2)=8
       GO TO 2051
              ENDIF 
                IF((XDIR.LE.0.8).AND.(XDIR.GT.0.6))THEN
                LPION(NNN,IRUN)=27
                EPION(NNN,IRUN)=Arho
                LB(I1)=9
                LB(I2)=6
       GO TO 2051
              ENDIF 
                IF(XDIR.GT.0.8)THEN
                LPION(NNN,IRUN)=27
                EPION(NNN,IRUN)=Arho
                LB(I1)=7
                LB(I2)=8
       GO TO 2051
              ENDIF 
               ENDIF
                IF(iabs(LB(I1)).EQ.2.AND.iabs(LB(I2)).EQ.2)THEN
                IF(XDIR.Le.0.2)then
                LPION(NNN,IRUN)=26
                EPION(NNN,IRUN)=Arho
              LB(I1)=6
              LB(I2)=7
       GO TO 2051
                ENDIF
                IF((XDIR.LE.0.4).AND.(XDIR.GT.0.2))THEN
                LPION(NNN,IRUN)=25
                EPION(NNN,IRUN)=Arho
                LB(I1)=6
                LB(I2)=9
       GO TO 2051
              ENDIF 
                IF((XDIR.GT.0.4).AND.(XDIR.LE.0.6))THEN
                LPION(NNN,IRUN)=27
                EPION(NNN,IRUN)=Arho
                LB(I1)=9
                LB(I2)=8
       GO TO 2051
              ENDIF 
                IF((XDIR.GT.0.6).AND.(XDIR.LE.0.8))THEN
                LPION(NNN,IRUN)=26
                EPION(NNN,IRUN)=Arho
                LB(I1)=7
                LB(I2)=7
       GO TO 2051
              ENDIF 
                IF(XDIR.GT.0.8)THEN
                LPION(NNN,IRUN)=25
                EPION(NNN,IRUN)=Arho
                LB(I1)=7
                LB(I2)=8
       GO TO 2051
              ENDIF 
              ENDIF
                IF(LB(I1)*LB(I2).EQ.2)THEN
                IF(XDIR.Le.0.17)then
                LPION(NNN,IRUN)=25
                EPION(NNN,IRUN)=Arho
              LB(I1)=6
              LB(I2)=9
       GO TO 2051
                ENDIF
                IF((XDIR.LE.0.34).AND.(XDIR.GT.0.17))THEN
                LPION(NNN,IRUN)=25
                EPION(NNN,IRUN)=Arho
                LB(I1)=7
                LB(I2)=9
       GO TO 2051
              ENDIF 
                IF((XDIR.GT.0.34).AND.(XDIR.LE.0.51))THEN
                LPION(NNN,IRUN)=27
                EPION(NNN,IRUN)=Arho
                LB(I1)=7
                LB(I2)=8
       GO TO 2051
              ENDIF 
                IF((XDIR.GT.0.51).AND.(XDIR.LE.0.68))THEN
                LPION(NNN,IRUN)=25
                EPION(NNN,IRUN)=Arho
                LB(I1)=8
                LB(I2)=8
       GO TO 2051
              ENDIF 
                IF((XDIR.GT.0.68).AND.(XDIR.LE.0.85))THEN
                LPION(NNN,IRUN)=26
                EPION(NNN,IRUN)=Arho
                LB(I1)=7
                LB(I2)=8
       GO TO 2051
              ENDIF 
                IF(XDIR.GT.0.85)THEN
                LPION(NNN,IRUN)=27
                EPION(NNN,IRUN)=Arho
                LB(I1)=7
                LB(I2)=7
              ENDIF 
                ENDIF
2051          E1CM    = SQRT (dm3**2 + PX3**2 + PY3**2 + PZ3**2)
              P1BETA  = PX3*BETAX + PY3*BETAY + PZ3*BETAZ
              TRANSF  = GAMMA * ( GAMMA * P1BETA / (GAMMA + 1) + E1CM )
              Pt1i1 = BETAX * TRANSF + PX3
              Pt2i1 = BETAY * TRANSF + PY3
              Pt3i1 = BETAZ * TRANSF + PZ3
             Eti1   = DM3
             if(ianti.eq.1 .and. lb(i1).ge.1 .and. lb(i2).ge.1)then
               lb(i1) = -lb(i1)
               lb(i2) = -lb(i2)
                if(LPION(NNN,IRUN) .eq. 25)then
                  LPION(NNN,IRUN)=27
                elseif(LPION(NNN,IRUN) .eq. 27)then
                  LPION(NNN,IRUN)=25
                endif
               endif
             lb1=lb(i1)
                E2CM    = SQRT (dm4**2 + PX4**2 + PY4**2 + PZ4**2)
                P2BETA  = PX4*BETAX+PY4*BETAY+PZ4*BETAZ
                TRANSF  = GAMMA * (GAMMA*P2BETA / (GAMMA + 1.) + E2CM)
                Pt1I2 = BETAX * TRANSF + PX4
                Pt2I2 = BETAY * TRANSF + PY4
                Pt3I2 = BETAZ * TRANSF + PZ4
              EtI2   = DM4
              lb2=lb(i2)
              p(1,i1)=pt1i1
              p(2,i1)=pt2i1
              p(3,i1)=pt3i1
              e(i1)=eti1
              lb(i1)=lb1
              p(1,i2)=pt1i2
              p(2,i2)=pt2i2
              p(3,i2)=pt3i2
              e(i2)=eti2
              lb(i2)=lb2
                PX1     = P(1,I1)
                PY1     = P(2,I1)
                PZ1     = P(3,I1)
              EM1       = E(I1)
                ID(I1)  = 2
                ID(I2)  = 2
                ID1     = ID(I1)
                IBLOCK=44
                EPCM=SQRT(EPION(NNN,IRUN)**2+PPX**2+PPY**2+PPZ**2)
                PPBETA=PPX*BETAX+PPY*BETAY+PPZ*BETAZ
                TRANSF=GAMMA*(GAMMA*PPBETA/(GAMMA+1.)+EPCM)
                PPION(1,NNN,IRUN)=BETAX*TRANSF+PPX
                PPION(2,NNN,IRUN)=BETAY*TRANSF+PPY
                PPION(3,NNN,IRUN)=BETAZ*TRANSF+PPZ
                dppion(nnn,irun)=dpertp(i1)*dpertp(i2)
                RPION(1,NNN,IRUN)=R(1,I1)
                RPION(2,NNN,IRUN)=R(2,I1)
                RPION(3,NNN,IRUN)=R(3,I1)
              go to 90005
308     CONTINUE
           NTRY1=0
126        CALL pprho(SRT,ISEED,PX3,PY3,PZ3,DM3,PX4,PY4,PZ4,DM4,
     &  PPX,PPY,PPZ,amrho,icou1)
       NTRY1=NTRY1+1
       if((icou1.lt.0).AND.(NTRY1.LE.20))GO TO 126
       CALL ROTATE(PX,PY,PZ,PX3,PY3,PZ3)
       CALL ROTATE(PX,PY,PZ,PX4,PY4,PZ4)
       CALL ROTATE(PX,PY,PZ,PPX,PPY,PPZ)
                NNN=NNN+1
              arho=amrho
              XDIR=RANART(NSEED)
                IF(LB(I1)*LB(I2).EQ.1)THEN
                IF(XDIR.Le.0.5)then
                LPION(NNN,IRUN)=26
                EPION(NNN,IRUN)=Arho
              LB(I1)=1
              LB(I2)=1
       GO TO 2052
                Else
                LPION(NNN,IRUN)=27
                EPION(NNN,IRUN)=Arho
                LB(I1)=1
                LB(I2)=2
       GO TO 2052
              ENDIF 
              endif
                IF(iabs(LB(I1)).EQ.2.AND.iabs(LB(I2)).EQ.2)THEN
                IF(XDIR.Le.0.5)then
                LPION(NNN,IRUN)=26
                EPION(NNN,IRUN)=Arho
              LB(I1)=2
              LB(I2)=2
       GO TO 2052
                Else
                LPION(NNN,IRUN)=25
                EPION(NNN,IRUN)=Arho
                LB(I1)=1
                LB(I2)=2
       GO TO 2052
              ENDIF 
              endif
                IF(LB(I1)*LB(I2).EQ.2)THEN
                IF(XDIR.Le.0.33)then
                LPION(NNN,IRUN)=26
                EPION(NNN,IRUN)=Arho
              LB(I1)=1
              LB(I2)=2
       GO TO 2052
                else IF((XDIR.LE.0.67).AND.(XDIR.GT.0.34))THEN
                LPION(NNN,IRUN)=25
                EPION(NNN,IRUN)=Arho
                LB(I1)=1
                LB(I2)=1
       GO TO 2052
              Else 
                LPION(NNN,IRUN)=27
                EPION(NNN,IRUN)=Arho
                LB(I1)=2
                LB(I2)=2
       GO TO 2052
              ENDIF 
              endif
2052          E1CM    = SQRT (dm3**2 + PX3**2 + PY3**2 + PZ3**2)
              P1BETA  = PX3*BETAX + PY3*BETAY + PZ3*BETAZ
              TRANSF  = GAMMA * ( GAMMA * P1BETA / (GAMMA + 1) + E1CM )
              Pt1i1 = BETAX * TRANSF + PX3
              Pt2i1 = BETAY * TRANSF + PY3
              Pt3i1 = BETAZ * TRANSF + PZ3
             Eti1   = DM3
              if(ianti.eq.1 .and. lb(i1).ge.1 .and. lb(i2).ge.1)then
               lb(i1) = -lb(i1)
               lb(i2) = -lb(i2)
                if(LPION(NNN,IRUN) .eq. 25)then
                  LPION(NNN,IRUN)=27
                elseif(LPION(NNN,IRUN) .eq. 27)then
                  LPION(NNN,IRUN)=25
                endif
               endif
             lb1=lb(i1)
                E2CM    = SQRT (dm4**2 + PX4**2 + PY4**2 + PZ4**2)
                P2BETA  = PX4*BETAX+PY4*BETAY+PZ4*BETAZ
                TRANSF  = GAMMA * (GAMMA*P2BETA / (GAMMA + 1.) + E2CM)
                Pt1I2 = BETAX * TRANSF + PX4
                Pt2I2 = BETAY * TRANSF + PY4
                Pt3I2 = BETAZ * TRANSF + PZ4
              EtI2   = DM4
              lb2=lb(i2)
              p(1,i1)=pt1i1
              p(2,i1)=pt2i1
              p(3,i1)=pt3i1
              e(i1)=eti1
              lb(i1)=lb1
              p(1,i2)=pt1i2
              p(2,i2)=pt2i2
              p(3,i2)=pt3i2
              e(i2)=eti2
              lb(i2)=lb2
                PX1     = P(1,I1)
                PY1     = P(2,I1)
                PZ1     = P(3,I1)
              EM1       = E(I1)
                ID(I1)  = 2
                ID(I2)  = 2
                ID1     = ID(I1)
                IBLOCK=45
                EPCM=SQRT(EPION(NNN,IRUN)**2+PPX**2+PPY**2+PPZ**2)
                PPBETA=PPX*BETAX+PPY*BETAY+PPZ*BETAZ
                TRANSF=GAMMA*(GAMMA*PPBETA/(GAMMA+1.)+EPCM)
                PPION(1,NNN,IRUN)=BETAX*TRANSF+PPX
                PPION(2,NNN,IRUN)=BETAY*TRANSF+PPY
                PPION(3,NNN,IRUN)=BETAZ*TRANSF+PPZ
                dppion(nnn,irun)=dpertp(i1)*dpertp(i2)
                RPION(1,NNN,IRUN)=R(1,I1)
                RPION(2,NNN,IRUN)=R(2,I1)
                RPION(3,NNN,IRUN)=R(3,I1)
              go to 90005
309     CONTINUE
           NTRY1=0
138        CALL ppomga(SRT,ISEED,PX3,PY3,PZ3,DM3,PX4,PY4,PZ4,DM4,
     &  PPX,PPY,PPZ,icou1)
       NTRY1=NTRY1+1
       if((icou1.lt.0).AND.(NTRY1.LE.20))GO TO 138
       CALL ROTATE(PX,PY,PZ,PX3,PY3,PZ3)
       CALL ROTATE(PX,PY,PZ,PX4,PY4,PZ4)
       CALL ROTATE(PX,PY,PZ,PPX,PPY,PPZ)
                NNN=NNN+1
              aomega=0.782
                IF(LB(I1)*LB(I2).EQ.1)THEN
                LPION(NNN,IRUN)=28
                EPION(NNN,IRUN)=Aomega
              LB(I1)=1
              LB(I2)=1
       GO TO 2053
                ENDIF
                IF(iabs(LB(I1)).EQ.2.AND.iabs(LB(I2)).EQ.2)THEN
                LPION(NNN,IRUN)=28
                EPION(NNN,IRUN)=Aomega
              LB(I1)=2
              LB(I2)=2
       GO TO 2053
                ENDIF
                IF(LB(I1)*LB(I2).EQ.2)THEN
                LPION(NNN,IRUN)=28
                EPION(NNN,IRUN)=Aomega
              LB(I1)=1
              LB(I2)=2
       GO TO 2053
                ENDIF
2053          E1CM    = SQRT (dm3**2 + PX3**2 + PY3**2 + PZ3**2)
              P1BETA  = PX3*BETAX + PY3*BETAY + PZ3*BETAZ
              TRANSF  = GAMMA * ( GAMMA * P1BETA / (GAMMA + 1) + E1CM )
              Pt1i1 = BETAX * TRANSF + PX3
              Pt2i1 = BETAY * TRANSF + PY3
              Pt3i1 = BETAZ * TRANSF + PZ3
             Eti1   = DM3
              if(ianti.eq.1 .and. lb(i1).ge.1 .and. lb(i2).ge.1)then
               lb(i1) = -lb(i1)
               lb(i2) = -lb(i2)
               endif
             lb1=lb(i1)
                E2CM    = SQRT (dm4**2 + PX4**2 + PY4**2 + PZ4**2)
                P2BETA  = PX4*BETAX+PY4*BETAY+PZ4*BETAZ
                TRANSF  = GAMMA * (GAMMA*P2BETA / (GAMMA + 1.) + E2CM)
                Pt1I2 = BETAX * TRANSF + PX4
                Pt2I2 = BETAY * TRANSF + PY4
                Pt3I2 = BETAZ * TRANSF + PZ4
              EtI2   = DM4
                lb2=lb(i2)
              p(1,i1)=pt1i1
              p(2,i1)=pt2i1
              p(3,i1)=pt3i1
              e(i1)=eti1
              lb(i1)=lb1
              p(1,i2)=pt1i2
              p(2,i2)=pt2i2
              p(3,i2)=pt3i2
              e(i2)=eti2
              lb(i2)=lb2
                PX1     = P(1,I1)
                PY1     = P(2,I1)
                PZ1     = P(3,I1)
              EM1       = E(I1)
                ID(I1)  = 2
                ID(I2)  = 2
                ID1     = ID(I1)
                IBLOCK=46
                EPCM=SQRT(EPION(NNN,IRUN)**2+PPX**2+PPY**2+PPZ**2)
                PPBETA=PPX*BETAX+PPY*BETAY+PPZ*BETAZ
                TRANSF=GAMMA*(GAMMA*PPBETA/(GAMMA+1.)+EPCM)
                PPION(1,NNN,IRUN)=BETAX*TRANSF+PPX
                PPION(2,NNN,IRUN)=BETAY*TRANSF+PPY
                PPION(3,NNN,IRUN)=BETAZ*TRANSF+PPZ
                dppion(nnn,irun)=dpertp(i1)*dpertp(i2)
                    RPION(1,NNN,IRUN)=R(1,I1)
                    RPION(2,NNN,IRUN)=R(2,I1)
                    RPION(3,NNN,IRUN)=R(3,I1)
              go to 90005
90005       continue
       RETURN
107     IF(PX .EQ. 0.0 .AND. PY .EQ. 0.0) THEN
        T2 = 0.0
      ELSE
        T2=ATAN2(PY,PX)
      END IF
      S1   = 1.0 - C1**2 
       IF(S1.LE.0)S1=0
       S1=SQRT(S1)
       scheck=1.0 - C2**2
       if(scheck.lt.0) then
          write(99,*) 'scheck3: ', scheck
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
      END
