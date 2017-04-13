      SUBROUTINE CRRD(PX,PY,PZ,SRT,I1,I2,
     & IBLOCK,xkaon0,xkaon,Xphi,xphin)
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,
     2  AP2=0.13957,AM0=1.232,PI=3.1415926,CUTOFF=1.8966,AVMASS=0.9383)
        PARAMETER     (AKA=0.498,ALA=1.1157,ASA=1.1974,APHI=1.02)
        parameter     (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
        COMMON /AA/ R(3,MAXSTR)
        COMMON /BB/ P(3,MAXSTR)
        COMMON /CC/ E(MAXSTR)
        COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
        common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      COMMON/RNDF77/NSEED
      SAVE   
       PX0=PX
       PY0=PY
       PZ0=PZ
       IBLOCK=1
       ianti=0
       if(lb(i1).lt.0 .or. lb(i2).lt.0) ianti=1
       x1=RANART(NSEED)
       if(xkaon0/(xkaon+Xphi).ge.x1)then
        IBLOCK=7
        if(ianti .eq. 1)iblock=-7
        NTAG=0
       KAONC=0
       IF(PNLKA(SRT)/(PNLKA(SRT)
     & +PNSKA(SRT)).GT.RANART(NSEED))KAONC=1
       IF(E(I1).LE.0.92)THEN
       LB(I1)=23
       E(I1)=AKA
              IF(KAONC.EQ.1)THEN
       LB(I2)=14
       E(I2)=ALA
              ELSE
        LB(I2) = 15 + int(3 * RANART(NSEED))
       E(I2)=ASA       
              ENDIF
         if(ianti .eq. 1)then
          lb(i1) = 21
          lb(i2) = -lb(i2)
         endif
       ELSE
       LB(I2)=23
       E(I2)=AKA
              IF(KAONC.EQ.1)THEN
       LB(I1)=14
       E(I1)=ALA
              ELSE
         LB(I1) = 15 + int(3 * RANART(NSEED))
       E(I1)=ASA       
              ENDIF
         if(ianti .eq. 1)then
          lb(i2) = 21
          lb(i1) = -lb(i1)
         endif
       ENDIF
        EM1=E(I1)
        EM2=E(I2)
       go to 50
       elseif(Xphi/(xkaon+Xphi).ge.x1)then
          iblock=222
         if(xphin/Xphi .ge. RANART(NSEED))then
          LB(I1)= 1+int(2*RANART(NSEED))
           E(I1)=AMN
         else
          LB(I1)= 6+int(4*RANART(NSEED))
           E(I1)=AM0
         endif
         if(ianti .eq. 1)lb(i1)=-lb(i1)
          LB(I2)= 29
           E(I2)=APHI
        EM1=E(I1)
        EM2=E(I2)
       go to 50
         else
       X2=RANART(NSEED)
       IBLOCK=81
       ntag=0
       if(lb(i1).eq.28.or.lb(i2).eq.28)go to 60
       if( ((lb(i1).eq.8.and.lb(i2).eq.27).
     &  or.(lb(i1).eq.27.and.lb(i2).eq.8))
     &        .OR. ((lb(i1).eq.-8.and.lb(i2).eq.25).
     &  or.(lb(i1).eq.25.and.lb(i2).eq.-8)) )then
              if(iabs(lb(i1)).eq.8)then
        ii = i1
       lb(i1)=1
       e(i1)=amn
       lb(i2)=5
       e(i2)=ap1
       go to 40
       else
        ii = i2
       lb(i2)=1
       e(i2)=amn
       lb(i1)=5
       e(i1)=ap1
       go to 40
              endif
       endif
       if((iabs(lb(i1)).eq.7.and.lb(i2).eq.26).
     &  or.(lb(i1).eq.26.and.iabs(lb(i2)).eq.7))then
              if(iabs(lb(i1)).eq.7)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=2
       e(i1)=amn
       lb(i2)=4
       e(i2)=ap1
       go to 40
       Else
       lb(i1)=1
       e(i1)=amn
       lb(i2)=3
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=2
       e(i2)=amn
       lb(i1)=4
       e(i1)=ap1
       go to 40
       Else
       lb(i2)=1
       e(i2)=amn
       lb(i1)=3
       e(i1)=ap1
       go to 40
       endif
              endif
       endif
       if((iabs(lb(i1)).eq.8.and.lb(i2).eq.26).
     &  or.(lb(i1).eq.26.and.iabs(lb(i2)).eq.8))then
              if(iabs(lb(i1)).eq.8)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=2
       e(i1)=amn
       lb(i2)=5
       e(i2)=ap1
       go to 40
       Else
       lb(i1)=1
       e(i1)=amn
       lb(i2)=4
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=2
       e(i2)=amn
       lb(i1)=5
       e(i1)=ap1
       go to 40
       Else
       lb(i2)=1
       e(i2)=amn
       lb(i1)=4
       e(i1)=ap1
       go to 40
       endif
              endif
       endif
       if((iabs(lb(i1)).eq.6.and.lb(i2).eq.26).
     &  or.(lb(i1).eq.26.and.iabs(lb(i2)).eq.6))then
              if(iabs(lb(i1)).eq.6)then
        ii = i1
       lb(i1)=2
       e(i1)=amn
       lb(i2)=3
       e(i2)=ap1
       go to 40
       else
        ii = i2
       lb(i2)=2
       e(i2)=amn
       lb(i1)=3
       e(i1)=ap1
       go to 40
       ENDIF
       endif
       if( ((lb(i1).eq.8.and.lb(i2).eq.25).
     &  or.(lb(i1).eq.25.and.lb(i2).eq.8))
     &        .OR. ((lb(i1).eq.-8.and.lb(i2).eq.27).
     &  or.(lb(i1).eq.27.and.lb(i2).eq.-8)) )then
              if(iabs(lb(i1)).eq.8)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=2
       e(i1)=amn
       lb(i2)=4
       e(i2)=ap1
       go to 40
       ELSE
       lb(i1)=1
       e(i1)=amn
       lb(i2)=3
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=2
       e(i2)=amn
       lb(i1)=4
       e(i1)=ap1
       go to 40
       ELSE
       lb(i2)=1
       e(i2)=amn
       lb(i1)=3
       e(i1)=ap1
       go to 40
       endif
              endif
       ENDIF
       if( ((lb(i1).eq.7.and.lb(i2).eq.27).
     &  or.(lb(i1).eq.27.and.lb(i2).eq.7))
     &       .OR.((lb(i1).eq.-7.and.lb(i2).eq.25).
     &  or.(lb(i1).eq.25.and.lb(i2).eq.-7)) )then
              if(iabs(lb(i1)).eq.7)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=2
       e(i1)=amn
       lb(i2)=5
       e(i2)=ap1
       go to 40
       else
       lb(i1)=1
       e(i1)=amn
       lb(i2)=4
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=2
       e(i2)=amn
       lb(i1)=5
       e(i1)=ap1
       go to 40
       Else
       lb(i2)=1
       e(i2)=amn
       lb(i1)=4
       e(i1)=ap1
       go to 40
       endif
              endif
       ENDIF
       if( ((lb(i1).eq.7.and.lb(i2).eq.25).
     &  or.(lb(i1).eq.25.and.lb(i2).eq.7))
     &       .OR.((lb(i1).eq.-7.and.lb(i2).eq.27).
     &  or.(lb(i1).eq.27.and.lb(i2).eq.-7)) )then
              if(iabs(lb(i1)).eq.7)then
        ii = i1
       lb(i1)=2
       e(i1)=amn
       lb(i2)=3
       e(i2)=ap1
       go to 40
       else
        ii = i2
       lb(i2)=2
       e(i2)=amn
       lb(i1)=3
       e(i1)=ap1
       go to 40
       ENDIF
       endif
       if( ((lb(i1).eq.6.and.lb(i2).eq.27).
     &  or.(lb(i1).eq.27.and.lb(i2).eq.6))
     &        .OR. ((lb(i1).eq.-6.and.lb(i2).eq.25).
     &  or.(lb(i1).eq.25.and.lb(i2).eq.-6)) )then
              if(iabs(lb(i1)).eq.6)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=2
       e(i1)=amn
       lb(i2)=4
       e(i2)=ap1
       go to 40
       else
       lb(i1)=1
       e(i1)=amn
       lb(i2)=3
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=2
       e(i2)=amn
       lb(i1)=4
       e(i1)=ap1
       go to 40
       Else
       lb(i2)=1
       e(i2)=amn
       lb(i1)=3
       e(i1)=ap1
       go to 40
       endif
              endif
       ENDIF
       if( ((lb(i1).eq.9.and.lb(i2).eq.25).
     &  or.(lb(i1).eq.25.and.lb(i2).eq.9))
     &        .OR.((lb(i1).eq.-9.and.lb(i2).eq.27).
     &  or.(lb(i1).eq.27.and.lb(i2).eq.-9)) )then
              if(iabs(lb(i1)).eq.9)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=2
       e(i1)=amn
       lb(i2)=5
       e(i2)=ap1
       go to 40
       else
       lb(i1)=1
       e(i1)=amn
       lb(i2)=4
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=2
       e(i2)=amn
       lb(i1)=5
       e(i1)=ap1
       go to 40
       Else
       lb(i2)=1
       e(i2)=amn
       lb(i1)=4
       e(i1)=ap1
       go to 40
       endif
              endif
       ENDIF
       if((iabs(lb(i1)).eq.9.and.lb(i2).eq.26).
     &  or.(lb(i1).eq.26.and.iabs(lb(i2)).eq.9))then
              if(iabs(lb(i1)).eq.9)then
        ii = i1
       lb(i1)=1
       e(i1)=amn
       lb(i2)=5
       e(i2)=ap1
       go to 40
       else
        ii = i2
       lb(i2)=1
       e(i2)=amn
       lb(i1)=5
       e(i1)=ap1
       go to 40
       ENDIF
       endif
       if( ((lb(i1).eq.11.and.lb(i2).eq.27).
     &  or.(lb(i1).eq.27.and.lb(i2).eq.11).
     &  or.(lb(i1).eq.13.and.lb(i2).eq.27).
     &  or.(lb(i1).eq.27.and.lb(i2).eq.13))
     &        .OR. ((lb(i1).eq.-11.and.lb(i2).eq.25).
     &  or.(lb(i1).eq.25.and.lb(i2).eq.-11).
     &  or.(lb(i1).eq.-13.and.lb(i2).eq.25).
     &  or.(lb(i1).eq.25.and.lb(i2).eq.-13)) )then
              if(iabs(lb(i1)).eq.11.or.iabs(lb(i1)).eq.13)then
        ii = i1
       lb(i1)=1
       e(i1)=amn
       lb(i2)=5
       e(i2)=ap1
       go to 40
       else
        ii = i2
       lb(i2)=1
       e(i2)=amn
       lb(i1)=5
       e(i1)=ap1
       go to 40
              endif
       endif
       if((iabs(lb(i1)).eq.10.and.lb(i2).eq.26).
     &  or.(lb(i1).eq.26.and.iabs(lb(i2)).eq.10).
     &  or.(lb(i1).eq.26.and.iabs(lb(i2)).eq.12).
     &  or.(lb(i2).eq.26.and.iabs(lb(i1)).eq.12))then
              if(iabs(lb(i1)).eq.10.or.iabs(lb(i1)).eq.12)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=2
       e(i1)=amn
       lb(i2)=4
       e(i2)=ap1
       go to 40
       Else
       lb(i1)=1
       e(i1)=amn
       lb(i2)=3
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=2
       e(i2)=amn
       lb(i1)=4
       e(i1)=ap1
       go to 40
       Else
       lb(i2)=1
       e(i2)=amn
       lb(i1)=3
       e(i1)=ap1
       go to 40
       endif
              endif
       endif
       if((iabs(lb(i1)).eq.11.and.lb(i2).eq.26).
     &  or.(lb(i1).eq.26.and.iabs(lb(i2)).eq.11).
     &  or.(lb(i1).eq.26.and.iabs(lb(i2)).eq.13).
     &  or.(lb(i2).eq.26.and.iabs(lb(i1)).eq.13))then
              if(iabs(lb(i1)).eq.11.or.iabs(lb(i1)).eq.13)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=2
       e(i1)=amn
       lb(i2)=5
       e(i2)=ap1
       go to 40
       Else
       lb(i1)=1
       e(i1)=amn
       lb(i2)=4
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=2
       e(i2)=amn
       lb(i1)=5
       e(i1)=ap1
       go to 40
       Else
       lb(i2)=1
       e(i2)=amn
       lb(i1)=4
       e(i1)=ap1
       go to 40
       endif
              endif
       endif
       if( ((lb(i1).eq.11.and.lb(i2).eq.25).
     &  or.(lb(i1).eq.25.and.lb(i2).eq.11).
     &  or.(lb(i1).eq.25.and.lb(i2).eq.13).
     &  or.(lb(i2).eq.25.and.lb(i1).eq.13))
     &        .OR.((lb(i1).eq.-11.and.lb(i2).eq.27).
     &  or.(lb(i1).eq.27.and.lb(i2).eq.-11).
     &  or.(lb(i1).eq.27.and.lb(i2).eq.-13).
     &  or.(lb(i2).eq.27.and.lb(i1).eq.-13)) )then
       if(iabs(lb(i1)).eq.11.or.iabs(lb(i1)).eq.13)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=2
       e(i1)=amn
       lb(i2)=4
       e(i2)=ap1
       go to 40
       ELSE
       lb(i1)=1
       e(i1)=amn
       lb(i2)=3
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=2
       e(i2)=amn
       lb(i1)=4
       e(i1)=ap1
       go to 40
       ELSE
       lb(i2)=1
       e(i2)=amn
       lb(i1)=3
       e(i1)=ap1
       go to 40
       endif
              endif
       ENDIF
       if( ((lb(i1).eq.10.and.lb(i2).eq.27).
     &  or.(lb(i1).eq.27.and.lb(i2).eq.10).
     &  or.(lb(i1).eq.12.and.lb(i2).eq.27).
     &  or.(lb(i1).eq.27.and.lb(i2).eq.12))
     &         .OR.((lb(i1).eq.-10.and.lb(i2).eq.25).
     &  or.(lb(i1).eq.25.and.lb(i2).eq.-10).
     &  or.(lb(i1).eq.-12.and.lb(i2).eq.25).
     &  or.(lb(i1).eq.25.and.lb(i2).eq.-12)) )then
       if(iabs(lb(i1)).eq.10.or.iabs(lb(i1)).eq.12)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=2
       e(i1)=amn
       lb(i2)=5
       e(i2)=ap1
       go to 40
       else
       lb(i1)=1
       e(i1)=amn
       lb(i2)=4
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=2
       e(i2)=amn
       lb(i1)=5
       e(i1)=ap1
       go to 40
       Else
       lb(i2)=1
       e(i2)=amn
       lb(i1)=4
       e(i1)=ap1
       go to 40
       endif
              endif
       ENDIF
       if( ((lb(i1).eq.10.and.lb(i2).eq.25).
     &  or.(lb(i1).eq.25.and.lb(i2).eq.10).
     &  or.(lb(i1).eq.25.and.lb(i2).eq.12).
     &  or.(lb(i1).eq.12.and.lb(i2).eq.25))
     &       .OR. ((lb(i1).eq.-10.and.lb(i2).eq.27).
     &  or.(lb(i1).eq.27.and.lb(i2).eq.-10).
     &  or.(lb(i1).eq.27.and.lb(i2).eq.-12).
     &  or.(lb(i1).eq.-12.and.lb(i2).eq.27)) )then
       if(iabs(lb(i1)).eq.10.or.iabs(lb(i1)).eq.12)then
        ii = i1
       lb(i1)=2
       e(i1)=amn
       lb(i2)=3
       e(i2)=ap1
       go to 40
       else
        ii = i2
       lb(i2)=2
       e(i2)=amn
       lb(i1)=3
       e(i1)=ap1
       go to 40
       ENDIF
       endif
60       IBLOCK=82
       if((iabs(lb(i1)).eq.7.and.lb(i2).eq.28).
     &  or.(lb(i1).eq.28.and.iabs(lb(i2)).eq.7))then
              if(iabs(lb(i1)).eq.7)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=2
       e(i1)=amn
       lb(i2)=4
       e(i2)=ap1
       go to 40
       Else
       lb(i1)=1
       e(i1)=amn
       lb(i2)=3
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=2
       e(i2)=amn
       lb(i1)=4
       e(i1)=ap1
       go to 40
       Else
       lb(i2)=1
       e(i2)=amn
       lb(i1)=3
       e(i1)=ap1
       go to 40
       endif
              endif
       endif
       if((iabs(lb(i1)).eq.8.and.lb(i2).eq.28).
     &  or.(lb(i1).eq.28.and.iabs(lb(i2)).eq.8))then
              if(iabs(lb(i1)).eq.8)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=2
       e(i1)=amn
       lb(i2)=5
       e(i2)=ap1
       go to 40
       Else
       lb(i1)=1
       e(i1)=amn
       lb(i2)=4
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=2
       e(i2)=amn
       lb(i1)=5
       e(i1)=ap1
       go to 40
       Else
       lb(i2)=1
       e(i2)=amn
       lb(i1)=4
       e(i1)=ap1
       go to 40
       endif
              endif
       endif
       if((iabs(lb(i1)).eq.6.and.lb(i2).eq.28).
     &  or.(lb(i1).eq.28.and.iabs(lb(i2)).eq.6))then
              if(iabs(lb(i1)).eq.6)then
        ii = i1
       lb(i1)=2
       e(i1)=amn
       lb(i2)=3
       e(i2)=ap1
       go to 40
       else
        ii = i2
       lb(i2)=2
       e(i2)=amn
       lb(i1)=3
       e(i1)=ap1
       go to 40
       ENDIF
       endif
       if((iabs(lb(i1)).eq.9.and.lb(i2).eq.28).
     &  or.(lb(i1).eq.28.and.iabs(lb(i2)).eq.9))then
              if(iabs(lb(i1)).eq.9)then
        ii = i1
       lb(i1)=1
       e(i1)=amn
       lb(i2)=5
       e(i2)=ap1
       go to 40
       else
        ii = i2
       lb(i2)=1
       e(i2)=amn
       lb(i1)=5
       e(i1)=ap1
       go to 40
       ENDIF
       endif
       if((iabs(lb(i1)).eq.10.and.lb(i2).eq.28).
     &  or.(lb(i1).eq.28.and.iabs(lb(i2)).eq.10).
     &  or.(lb(i1).eq.28.and.iabs(lb(i2)).eq.12).
     &  or.(lb(i2).eq.28.and.iabs(lb(i1)).eq.12))then
              if(iabs(lb(i1)).eq.10.or.iabs(lb(i1)).eq.12)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=2
       e(i1)=amn
       lb(i2)=4
       e(i2)=ap1
       go to 40
       Else
       lb(i1)=1
       e(i1)=amn
       lb(i2)=3
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=2
       e(i2)=amn
       lb(i1)=4
       e(i1)=ap1
       go to 40
       Else
       lb(i2)=1
       e(i2)=amn
       lb(i1)=3
       e(i1)=ap1
       go to 40
       endif
              endif
       endif
       if((iabs(lb(i1)).eq.11.and.lb(i2).eq.28).
     &  or.(lb(i1).eq.28.and.iabs(lb(i2)).eq.11).
     &  or.(lb(i1).eq.28.and.iabs(lb(i2)).eq.13).
     &  or.(lb(i2).eq.28.and.iabs(lb(i1)).eq.13))then
              if(iabs(lb(i1)).eq.11.or.iabs(lb(i1)).eq.13)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=2
       e(i1)=amn
       lb(i2)=5
       e(i2)=ap1
       go to 40
       Else
       lb(i1)=1
       e(i1)=amn
       lb(i2)=4
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=2
       e(i2)=amn
       lb(i1)=5
       e(i1)=ap1
       go to 40
       Else
       lb(i2)=1
       e(i2)=amn
       lb(i1)=4
       e(i1)=ap1
       go to 40
       endif
              endif
       endif
40       em1=e(i1)
       em2=e(i2)
       if(ianti.eq.1 .and. lb(i1).ge.1 .and. lb(i2).ge.1)then
         lb(ii) = -lb(ii)
           jj = i2
          if(ii .eq. i2)jj = i1
          if(lb(jj).eq.3)then
           lb(jj) = 5
          elseif(lb(jj).eq.5)then
           lb(jj) = 3
          endif
         endif
       endif
50          PR2   = (SRT**2 - EM1**2 - EM2**2)**2
     1                - 4.0 * (EM1*EM2)**2
          IF(PR2.LE.0.)PR2=1.E-09
          PR=SQRT(PR2)/(2.*SRT)
          xptr=0.33*pr
         cc1=ptr(xptr,iseed)
         scheck=pr**2-cc1**2
         if(scheck.lt.0) then
            write(99,*) 'scheck39: ', scheck
            scheck=0.
         endif
         c1=sqrt(scheck)/pr
          T1   = 2.0 * PI * RANART(NSEED)
      S1   = SQRT( 1.0 - C1**2 )
      CT1  = COS(T1)
      ST1  = SIN(T1)
      PZ   = PR * C1
      PX   = PR * S1*CT1 
      PY   = PR * S1*ST1 
       CALL ROTATE(PX0,PY0,PZ0,PX,PY,PZ)
      RETURN
      END
