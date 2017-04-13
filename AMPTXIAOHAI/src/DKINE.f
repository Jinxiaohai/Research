        SUBROUTINE DKINE(IRUN,I,NNN,NLAB,ISEED,wid,nt)
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
      COMMON/RNDF77/NSEED
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
        EXTERNAL IARFLV, INVFLV
      SAVE   
        PX=P(1,I)
        PY=P(2,I)
        PZ=P(3,I)
        RX=R(1,I)
        RY=R(2,I)
        RZ=R(3,I)
        DM=E(I)
        EDELTA=SQRT(DM**2+PX**2+PY**2+PZ**2)
        PM=EPION(NNN,IRUN)
        AM=AMP
        IF(NLAB.EQ.2)AM=AMN
        Q2=((DM**2-AM**2+PM**2)/(2.*DM))**2-PM**2
        IF(Q2.LE.0.)Q2=1.e-09
        Q=SQRT(Q2)
11      QX=1.-2.*RANART(NSEED)
        QY=1.-2.*RANART(NSEED)
        QZ=1.-2.*RANART(NSEED)
        QS=QX**2+QY**2+QZ**2
        IF(QS.GT.1.) GO TO 11
        PXP=Q*QX/SQRT(QS)
        PYP=Q*QY/SQRT(QS)
        PZP=Q*QZ/SQRT(QS)
        EP=SQRT(Q**2+PM**2)
        PXN=-PXP
        PYN=-PYP
        PZN=-PZP
        EN=SQRT(Q**2+AM**2)
        GD=EDELTA/DM
        FGD=GD/(1.+GD)
        BDX=PX/EDELTA
        BDY=PY/EDELTA
        BDZ=PZ/EDELTA
        BPP=BDX*PXP+BDY*PYP+BDZ*PZP
        BPN=BDX*PXN+BDY*PYN+BDZ*PZN
        P(1,I)=PXN+BDX*GD*(FGD*BPN+EN)
        P(2,I)=PYN+BDY*GD*(FGD*BPN+EN)
        P(3,I)=PZN+BDZ*GD*(FGD*BPN+EN)
        E(I)=AM
        PPION(1,NNN,IRUN)=PXP+BDX*GD*(FGD*BPP+EP)
        PPION(2,NNN,IRUN)=PYP+BDY*GD*(FGD*BPP+EP)
        PPION(3,NNN,IRUN)=PZP+BDZ*GD*(FGD*BPP+EP)
        dppion(NNN,IRUN)=dpertp(I)
        RPION(1,NNN,IRUN)=R(1,I)
        RPION(2,NNN,IRUN)=R(2,I)
        RPION(3,NNN,IRUN)=R(3,I)
        devio=SQRT(EPION(NNN,IRUN)**2+PPION(1,NNN,IRUN)**2
     1       +PPION(2,NNN,IRUN)**2+PPION(3,NNN,IRUN)**2)
     2       +SQRT(E(I)**2+P(1,I)**2+P(2,I)**2+P(3,I)**2)-e1
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
        endif
 200    format(a30,2(1x,e10.4))
 210    format(i6,5(1x,f8.3))
 220    format(a2,i5,5(1x,f8.3))
        RETURN
        END
