      SUBROUTINE XND(px,py,pz,srt,I1,I2,xinel,
     &               sigk,xsk1,xsk2,xsk3,xsk4,xsk5)
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
       xinel=0.
       sigk=0
       xsk1=0
       xsk2=0
       xsk3=0
       xsk4=0
       xsk5=0
        EM1=E(I1)
        EM2=E(I2)
      PR  = SQRT( PX**2 + PY**2 + PZ**2 )
        IF (SRT .LT. 2.04) RETURN
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
222       SIGK=XSK1+XSK2+XSK3+XSK4
        XSK1 = 2.0 * XSK1
        XSK2 = 2.0 * XSK2
        XSK3 = 2.0 * XSK3
        XSK4 = 2.0 * XSK4
        SIGK = 2.0 * SIGK + xsk5
       if(((iabs(lb(i1)).eq.2).and.(iabs(lb(i2)).eq.6)).OR. 
     &         ((iabs(lb(i2)).eq.2).and.(iabs(lb(i1)).eq.6)).OR.
     &         ((iabs(lb(i1)).eq.1).and.(iabs(lb(i2)).eq.9)).OR.
     &         ((iabs(lb(i2)).eq.1).and.(iabs(lb(i1)).eq.9)))THEN
       xinel=sigk
       return
       ENDIF
        IF(LB(I1)*LB(I2).EQ.18.AND.
     &    (iabs(LB(I1)).EQ.2.OR.iabs(LB(I2)).EQ.2))then
        SIGND=SIGMA(SRT,1,1,0)+0.5*SIGMA(SRT,1,1,1)
        SIGDN=0.25*SIGND*RENOM
        xinel=SIGDN+X1440+X1535+SIGK
       RETURN
       endif
        IF(LB(I1)*LB(I2).EQ.6.AND.
     &    (iabs(LB(I1)).EQ.1.OR.iabs(LB(I2)).EQ.1))THEN
        SIGND=SIGMA(SRT,1,1,0)+0.5*SIGMA(SRT,1,1,1)
        SIGDN=0.25*SIGND*RENOM
        xinel=SIGDN+X1440+X1535+SIGK
       RETURN
       endif
        IF(LB(I1)*LB(I2).EQ.8.AND.
     &    (iabs(LB(I1)).EQ.1.OR.iabs(LB(I2)).EQ.1))THEN
        SIGND=1.5*SIGMA(SRT,1,1,1)
        SIGDN=0.25*SIGND*RENOM
        xinel=SIGDN+x1440+x1535+SIGK
       RETURN
       endif
        IF(LB(I1)*LB(I2).EQ.14.AND.
     &   (iabs(LB(I1)).EQ.2.AND.iabs(LB(I2)).EQ.2))THEN
        SIGND=1.5*SIGMA(SRT,1,1,1)
        SIGDN=0.25*SIGND*RENOM
        xinel=SIGDN+x1440+x1535+SIGK
       RETURN
       endif
        IF(LB(I1)*LB(I2).EQ.16.AND.
     &     (iabs(LB(I1)).EQ.2.OR.iabs(LB(I2)).EQ.2))THEN
        SIGND=0.5*SIGMA(SRT,1,1,1)+0.25*SIGMA(SRT,1,1,0)
        SIGDN=0.5*SIGND*RENOM
        xinel=SIGDN+2.*x1440+2.*x1535+SIGK
       RETURN
       endif
        IF(LB(I1)*LB(I2).EQ.7)THEN
        SIGND=0.5*SIGMA(SRT,1,1,1)+0.25*SIGMA(SRT,1,1,0)
        SIGDN=0.5*SIGND*RENOM
        xinel=SIGDN+2.*x1440+2.*x1535+SIGK
       RETURN
       endif
        IF(LB(I1)*LB(I2).EQ.10.AND.
     &   (iabs(LB(I1)).EQ.1.OR.iabs(LB(I2)).EQ.1))then
        SIGND=(3./4.)*SIGMA(SRT,2,0,1)
        SIGDN=SIGND*RENOMN
        xinel=SIGDN+X1535+SIGK
       RETURN
       endif
        IF(LB(I1)*LB(I2).EQ.22.AND.
     &   (iabs(LB(I1)).EQ.2.OR.iabs(LB(I2)).EQ.2))then
        SIGND=(3./4.)*SIGMA(SRT,2,0,1)
        SIGDN=SIGND*RENOMN
        xinel=SIGDN+X1535+SIGK
       RETURN
       endif
        IF((iabs(LB(I1)).EQ.12).OR.(iabs(LB(I1)).EQ.13).OR.
     1  (iabs(LB(I2)).EQ.12).OR.(iabs(LB(I2)).EQ.13))THEN
        SIGND=X1535
        SIGDN=SIGND*RENOM1
        xinel=SIGDN+SIGK
       RETURN
       endif
        RETURN
       end
