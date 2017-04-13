        SUBROUTINE DKINE2(IRUN,I,NNN,NLAB,ISEED,wid,nt)
        PARAMETER (hbarc=0.19733)
        PARAMETER (MAXSTR=150001,MAXR=1,
     1  AMN=0.939457,AMP=0.93828,ETAM=0.5475,
     2  AP1=0.13496,AP2=0.13957,AM0=1.232,PI=3.1415926)
        COMMON /AA/ R(3,MAXSTR)
        COMMON /BB/ P(3,MAXSTR)
        COMMON /CC/ E(MAXSTR)
        COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
        COMMON   /RUN/NUM
        COMMON   /PA/RPION(3,MAXSTR,MAXR)
        COMMON   /PB/PPION(3,MAXSTR,MAXR)
        COMMON   /PC/EPION(MAXSTR,MAXR)
        COMMON   /PD/LPION(MAXSTR,MAXR)
      common/leadng/lb1,px1,py1,pz1,em1,e1,xfnl,yfnl,zfnl,tfnl,
     1 px1n,py1n,pz1n,dp1n
        COMMON/tdecay/tfdcy(MAXSTR),tfdpi(MAXSTR,MAXR),tft(MAXSTR)
        COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &       IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
        EXTERNAL IARFLV, INVFLV
      COMMON/RNDF77/NSEED
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
      SAVE   
        PX=P(1,I)
        PY=P(2,I)
        PZ=P(3,I)
        RX=R(1,I)
        RY=R(2,I)
        RZ=R(3,I)
        DM=E(I)
        EDELTA=SQRT(DM**2+PX**2+PY**2+PZ**2)
        PM1=EPION(NNN,IRUN)
        PM2=EPION(NNN+1,IRUN)
        AM=AMN
       IF(NLAB.EQ.1)AM=AMP
       PMAX2=(DM**2-(AM+PM1+PM2)**2)*(DM**2-(AM-PM1-PM2)**2)/4/DM**2
       scheck=PMAX2
       if(scheck.lt.0) then
          write(99,*) 'scheck15: ', scheck
          scheck=0.
       endif
       PMAX=SQRT(scheck)
       CSS=1.-2.*RANART(NSEED)
       SSS=SQRT(1-CSS**2)
       FAI=2*PI*RANART(NSEED)
       PX0=PMAX*SSS*COS(FAI)
       PY0=PMAX*SSS*SIN(FAI)
       PZ0=PMAX*CSS
       EP0=SQRT(PX0**2+PY0**2+PZ0**2+AM**2)
       BETAX=-PX0/(DM-EP0)
       BETAY=-PY0/(DM-EP0)
       BETAZ=-PZ0/(DM-EP0)
       scheck=1-BETAX**2-BETAY**2-BETAZ**2
       if(scheck.le.0) then
          write(99,*) 'scheck16: ', scheck
          stop
       endif
       GD1=1./SQRT(scheck)
       FGD1=GD1/(1+GD1)
        Q2=((DM-EP0)/(2.*GD1))**2-PM1**2
        IF(Q2.LE.0.)Q2=1.E-09
        Q=SQRT(Q2)
11      QX=1.-2.*RANART(NSEED)
        QY=1.-2.*RANART(NSEED)
        QZ=1.-2.*RANART(NSEED)
        QS=QX**2+QY**2+QZ**2
        IF(QS.GT.1.) GO TO 11
        PXP=Q*QX/SQRT(QS)
        PYP=Q*QY/SQRT(QS)
        PZP=Q*QZ/SQRT(QS)
        EP=SQRT(Q**2+PM1**2)
        PXN=-PXP
        PYN=-PYP
        PZN=-PZP
        EN=SQRT(Q**2+PM2**2)
        BPP1=BETAX*PXP+BETAY*PYP+BETAZ*PZP
        BPN1=BETAX*PXN+BETAY*PYN+BETAZ*PZN
        P1M=PXN+BETAX*GD1*(FGD1*BPN1+EN)
        P2M=PYN+BETAY*GD1*(FGD1*BPN1+EN)
        P3M=PZN+BETAZ*GD1*(FGD1*BPN1+EN)
       EPN=SQRT(P1M**2+P2M**2+P3M**2+PM2**2)
        P1P=PXP+BETAX*GD1*(FGD1*BPP1+EP)
        P2P=PYP+BETAY*GD1*(FGD1*BPP1+EP)
        P3P=PZP+BETAZ*GD1*(FGD1*BPP1+EP)
       EPP=SQRT(P1P**2+P2P**2+P3P**2+PM1**2)
        GD=EDELTA/DM
        FGD=GD/(1.+GD)
        BDX=PX/EDELTA
        BDY=PY/EDELTA
        BDZ=PZ/EDELTA
       BP0=BDX*PX0+BDY*PY0+BDZ*PZ0
        BPP=BDX*P1P+BDY*P2P+BDZ*P3P
        BPN=BDX*P1M+BDY*P2M+BDZ*P3M
        P(1,I)=PX0+BDX*GD*(FGD*BP0+EP0)
        P(2,I)=PY0+BDY*GD*(FGD*BP0+EP0)
        P(3,I)=PZ0+BDZ*GD*(FGD*BP0+EP0)
       E(I)=am
       ID(I)=0
       enucl=sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+e(i)**2)
        PPION(1,NNN,IRUN)=P1P+BDX*GD*(FGD*BPP+EPP)
        PPION(2,NNN,IRUN)=P2P+BDY*GD*(FGD*BPP+EPP)
        PPION(3,NNN,IRUN)=P3P+BDZ*GD*(FGD*BPP+EPP)
       epion1=sqrt(ppion(1,nnn,irun)**2
     &  +ppion(2,nnn,irun)**2+ppion(3,nnn,irun)**2
     &  +epion(nnn,irun)**2)
        RPION(1,NNN,IRUN)=R(1,I)
        RPION(2,NNN,IRUN)=R(2,I)
        RPION(3,NNN,IRUN)=R(3,I)
        PPION(1,NNN+1,IRUN)=P1M+BDX*GD*(FGD*BPN+EPN)
        PPION(2,NNN+1,IRUN)=P2M+BDY*GD*(FGD*BPN+EPN)
        PPION(3,NNN+1,IRUN)=P3M+BDZ*GD*(FGD*BPN+EPN)
        dppion(NNN,IRUN)=dpertp(I)
        dppion(NNN+1,IRUN)=dpertp(I)
       epion2=sqrt(ppion(1,nnn+1,irun)**2
     &  +ppion(2,nnn+1,irun)**2+ppion(3,nnn+1,irun)**2
     &  +epion(nnn+1,irun)**2)
        RPION(1,NNN+1,IRUN)=R(1,I)
        RPION(2,NNN+1,IRUN)=R(2,I)
        RPION(3,NNN+1,IRUN)=R(3,I)
        devio=SQRT(EPION(NNN,IRUN)**2+PPION(1,NNN,IRUN)**2
     1       +PPION(2,NNN,IRUN)**2+PPION(3,NNN,IRUN)**2)
     2       +SQRT(E(I)**2+P(1,I)**2+P(2,I)**2+P(3,I)**2)
     3       +SQRT(EPION(NNN+1,IRUN)**2+PPION(1,NNN+1,IRUN)**2
     4       +PPION(2,NNN+1,IRUN)**2+PPION(3,NNN+1,IRUN)**2)-e1
        if(nt.eq.ntmax) then
           tau0=hbarc/wid
           taudcy=tau0*(-1.)*alog(1.-RANART(NSEED))
           taudcy=taudcy*e1/em1
           tfnl=tfnl+taudcy
           xfnl=xfnl+px1/e1*taudcy
           yfnl=yfnl+py1/e1*taudcy
           zfnl=zfnl+pz1/e1*taudcy
           R(1,I)=xfnl
           R(2,I)=yfnl
           R(3,I)=zfnl
           tfdcy(I)=tfnl
           RPION(1,NNN,IRUN)=xfnl
           RPION(2,NNN,IRUN)=yfnl
           RPION(3,NNN,IRUN)=zfnl
           tfdpi(NNN,IRUN)=tfnl
           RPION(1,NNN+1,IRUN)=xfnl
           RPION(2,NNN+1,IRUN)=yfnl
           RPION(3,NNN+1,IRUN)=zfnl
           tfdpi(NNN+1,IRUN)=tfnl
        endif
 200    format(a30,2(1x,e10.4))
 210    format(i6,5(1x,f8.3))
 220    format(a2,i5,5(1x,f8.3))
        RETURN
        END
