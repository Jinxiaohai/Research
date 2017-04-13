      SUBROUTINE XDDIN(PX,PY,PZ,SRT,I1,I2,
     &XINEL,SIGK,XSK1,XSK2,XSK3,XSK4,XSK5)
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,AKA=0.498,APHI=1.020,
     2  AP2=0.13957,AM0=1.232,PI=3.1415926,CUTOFF=1.8966,AVMASS=0.9383)
        parameter     (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
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
      SAVE   
       XINEL=0
       SIGK=0
       XSK1=0
       XSK2=0
       XSK3=0
       XSK4=0
       XSK5=0
        EM1=E(I1)
        EM2=E(I2)
      PR  = SQRT( PX**2 + PY**2 + PZ**2 )
       call N1535(iabs(lb(i1)),iabs(lb(i2)),srt,X1535)
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
222       SIGK=XSK1+XSK2+XSK3+XSK4
        XSK1 = 2.0 * XSK1
        XSK2 = 2.0 * XSK2
        XSK3 = 2.0 * XSK3
        XSK4 = 2.0 * XSK4
        SIGK = 2.0 * SIGK + xsk5
        IDD=iabs(LB(I1)*LB(I2))
       s2d=reab2d(i1,i2,srt)
        S2D = 0.
       if(((iabs(lb(i1)).ge.12).and.(iabs(lb(i2)).ge.12)).OR.
     &       ((iabs(lb(i1)).ge.12).and.(iabs(lb(i2)).ge.6)).OR.
     &       ((iabs(lb(i2)).ge.12).and.(iabs(lb(i1)).ge.6)))THEN
       XINEL=sigk+s2d
       RETURN
       ENDIF
        IF((IDD.EQ.63).OR.(IDD.EQ.64).OR.(IDD.EQ.48).
     1  OR.(IDD.EQ.49).OR.(IDD.EQ.11*11).OR.(IDD.EQ.10*10).
     2  OR.(IDD.EQ.88).OR.(IDD.EQ.66).
     3  OR.(IDD.EQ.90).OR.(IDD.EQ.70))THEN
        XINEL=X1535+SIGK+s2d
       RETURN
        ENDIF
       IF((IDD.EQ.110).OR.(IDD.EQ.77).OR.(IDD.EQ.80))THEN
       XINEL=X1535+SIGK+s2d
       RETURN
       ENDIF       
       IF((IDD.EQ.54).OR.(IDD.EQ.56))THEN
        SIG2=(3./4.)*SIGMA(SRT,2,0,1)
        XINEL=2.*(SIG2+X1535)+SIGK+s2d
       RETURN
       ENDIF
       RETURN
       END
