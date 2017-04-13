      SUBROUTINE CRPN(PX,PY,PZ,SRT,I1,I2,
     & IBLOCK,xkaon0,xkaon,Xphi,xphin)
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,APHI=1.020,
     2  AP2=0.13957,AM0=1.232,PI=3.1415926,CUTOFF=1.8966,AVMASS=0.9383)
        PARAMETER      (AKA=0.498,ALA=1.1157,ASA=1.1974)
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
      iblock=1
      x1=RANART(NSEED)
      ianti=0
      if(lb(i1).lt.0 .or. lb(i2).lt.0) ianti=1
      if(xkaon0/(xkaon+Xphi).ge.x1)then
        IBLOCK=7
        if(ianti .eq. 1)iblock=-7
        NTAG=0
       KAONC=0
       IF(PNLKA(SRT)/(PNLKA(SRT)
     &       +PNSKA(SRT)).GT.RANART(NSEED))KAONC=1
       IF(E(I1).LE.0.2)THEN
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
       IF(RANART(NSEED).LE.TWOPI(SRT)/
     &  (TWOPI(SRT)+THREPI(SRT)+FOURPI(SRT)))THEN
       iblock=77
       ELSE 
        IF(THREPI(SRT)/(THREPI(SRT)+FOURPI(SRT)).
     &  GT.RANART(NSEED))THEN
       IBLOCK=78
       ELSE
       IBLOCK=79
       ENDIF
       endif
       ntag=0
       X2=RANART(NSEED)
       if(iblock.eq.77)then
       dmax=srt-ap1-0.02
       dm=rmass(dmax,iseed)
       if( ((lb(i1).eq.1.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.1))
     &       .OR. ((lb(i1).eq.-1.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.-1)) )then
              if(iabs(lb(i1)).eq.1)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=8
       e(i1)=dm
       lb(i2)=5
       e(i2)=ap1
       go to 40
       ELSE
       lb(i1)=9
       e(i1)=dm
       lb(i2)=4
        ipi = 4
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=8
       e(i2)=dm
       lb(i1)=5
       e(i1)=ap1
       go to 40
       ELSE
       lb(i2)=9
       e(i2)=dm
       lb(i1)=4
       e(i1)=ap1
       go to 40
       endif
              endif
       endif
       if( ((lb(i1).eq.1.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.1))
     &        .OR. ((lb(i1).eq.-1.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.-1)) )then
              if(iabs(lb(i1)).eq.1)then
        ii = i1
       IF(X2.LE.0.33)THEN
       lb(i1)=6
       e(i1)=dm
       lb(i2)=5
       e(i2)=ap1
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i1)=7
       e(i1)=dm
       lb(i2)=4
       e(i2)=ap1
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i1)=8
       e(i1)=dm
       lb(i2)=3
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.33)THEN
       lb(i2)=6
       e(i2)=dm
       lb(i1)=5
       e(i1)=ap1
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i2)=7
       e(i2)=dm
       lb(i1)=4
       e(i1)=ap1
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i2)=8
       e(i2)=dm
       lb(i1)=3
       e(i1)=ap1
       go to 40
       endif
              endif
       endif
       if( ((lb(i1).eq.2.and.lb(i2).eq.5).
     &   or.(lb(i1).eq.5.and.lb(i2).eq.2))
     & .OR. ((lb(i1).eq.-2.and.lb(i2).eq.3).
     &   or.(lb(i1).eq.3.and.lb(i2).eq.-2)) )then
              if(iabs(lb(i1)).eq.2)then
        ii = i1
       IF(X2.LE.0.33)THEN
       lb(i1)=8
       e(i1)=dm
       lb(i2)=4
       e(i2)=ap1
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i1)=7
       e(i1)=dm
       lb(i2)=5
       e(i2)=ap1
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i1)=9
       e(i1)=dm
       lb(i2)=3
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.33)THEN
       lb(i2)=8
       e(i2)=dm
       lb(i1)=4
       e(i1)=ap1
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i2)=7
       e(i2)=dm
       lb(i1)=5
       e(i1)=ap1
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i2)=9
       e(i2)=dm
       lb(i1)=3
       e(i1)=ap1
       go to 40
       endif
              endif
       endif
       if((iabs(lb(i1)).eq.1.and.lb(i2).eq.4).
     &  or.(lb(i1).eq.4.and.iabs(lb(i2)).eq.1))then
              if(iabs(lb(i1)).eq.1)then
        ii = i1
       IF(X2.LE.0.33)THEN
       lb(i1)=8
       e(i1)=dm
       lb(i2)=4
       e(i2)=ap1
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i1)=7
       e(i1)=dm
       lb(i2)=5
       e(i2)=ap1
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i1)=9
       e(i1)=dm
       lb(i2)=3
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.33)THEN
       lb(i2)=8
       e(i2)=dm
       lb(i1)=4
       e(i1)=ap1
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i2)=7
       e(i2)=dm
       lb(i1)=5
       e(i1)=ap1
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i2)=9
       e(i2)=dm
       lb(i1)=3
       e(i1)=ap1
       go to 40
       endif
              endif
       endif 
       if( ((lb(i1).eq.2.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.2))
     &         .OR. ((lb(i1).eq.-2.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.-2)) )then
              if(iabs(lb(i1)).eq.2)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=6
       e(i1)=dm
       lb(i2)=4
       e(i2)=ap1
       go to 40
       ELSE
       lb(i1)=7
       e(i1)=dm
       lb(i2)=3
       e(i2)=ap1
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=6
       e(i2)=dm
       lb(i1)=4
       e(i1)=ap1
       go to 40
       ELSE
       lb(i2)=7
       e(i2)=dm
       lb(i1)=3
       e(i1)=ap1
       go to 40
       endif
              endif
       ENDIF
       if((iabs(lb(i1)).eq.2.and.lb(i2).eq.4).
     &  or.(lb(i1).eq.4.and.iabs(lb(i2)).eq.2))then
              if(iabs(lb(i1)).eq.2)then
        ii = i1
       IF(X2.LE.0.33)THEN
       lb(i1)=7
       e(i1)=dm
       lb(i2)=4
       e(i2)=ap1
       go to 40
       Endif
       IF(X2.LE.0.67.AND.X2.GT.0.33)THEN       
       lb(i1)=6
       e(i1)=dm
       lb(i2)=5
       e(i2)=ap1
       go to 40
       endif
       IF(X2.GT.0.67)THEN
       LB(I1)=8
       E(I1)=DM
       LB(I2)=3
       E(I2)=AP1
       GO TO 40
       ENDIF
              else
        ii = i2
       IF(X2.LE.0.33)THEN
       lb(i2)=7
       e(i2)=dm
       lb(i1)=4
       e(i1)=ap1
       go to 40
       ENDIF
       IF(X2.LE.0.67.AND.X2.GT.0.33)THEN       
       lb(i2)=6
       e(i2)=dm
       lb(i1)=5
       e(i1)=ap1
       go to 40
       endif
       IF(X2.GT.0.67)THEN
       LB(I2)=8
       E(I2)=DM
       LB(I1)=3
       E(I1)=AP1
       GO TO 40
       ENDIF
              endif
       endif
                     ENDIF
       if(iblock.eq.78)then
       call Rmasdd(srt,1.232,0.77,1.08,
     &  0.28,ISEED,4,dm,ameson)
       arho=AMESON
       if( ((lb(i1).eq.1.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.1))
     &        .OR. ((lb(i1).eq.-1.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.-1)) )then
              if(iabs(lb(i1)).eq.1)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=8
       e(i1)=dm
       lb(i2)=27
       e(i2)=arho
       go to 40
       ELSE
       lb(i1)=9
       e(i1)=dm
       lb(i2)=26
       e(i2)=arho
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=8
       e(i2)=dm
       lb(i1)=27
       e(i1)=arho
       go to 40
       ELSE
       lb(i2)=9
       e(i2)=dm
       lb(i1)=26
       e(i1)=arho
       go to 40
       endif
              endif
       endif
       if( ((lb(i1).eq.1.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.1))
     &        .OR. ((lb(i1).eq.-1.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.-1)) )then
              if(iabs(lb(i1)).eq.1)then
        ii = i1
       IF(X2.LE.0.33)THEN
       lb(i1)=6
       e(i1)=dm
       lb(i2)=27
       e(i2)=arho
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i1)=7
       e(i1)=dm
       lb(i2)=26
       e(i2)=arho
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i1)=8
       e(i1)=dm
       lb(i2)=25
       e(i2)=arho
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.33)THEN
       lb(i2)=6
       e(i2)=dm
       lb(i1)=27
       e(i1)=arho
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i2)=7
       e(i2)=dm
       lb(i1)=26
       e(i1)=arho
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i2)=8
       e(i2)=dm
       lb(i1)=25
       e(i1)=arho
       go to 40
       endif
              endif
       endif
       if( ((lb(i1).eq.2.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.2))
     &       .OR.((lb(i1).eq.-2.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.-2)) )then
              if(iabs(lb(i1)).eq.2)then
        ii = i1
       IF(X2.LE.0.33)THEN
       lb(i1)=8
       e(i1)=dm
       lb(i2)=26
       e(i2)=arho
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i1)=7
       e(i1)=dm
       lb(i2)=27
       e(i2)=arho
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i1)=9
       e(i1)=dm
       lb(i2)=25
       e(i2)=arho
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.33)THEN
       lb(i2)=8
       e(i2)=dm
       lb(i1)=26
       e(i1)=arho
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i2)=7
       e(i2)=dm
       lb(i1)=27
       e(i1)=arho
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i2)=9
       e(i2)=dm
       lb(i1)=25
       e(i1)=arho
       go to 40
       endif
              endif
       endif
       if((iabs(lb(i1)).eq.1.and.lb(i2).eq.4).
     &  or.(lb(i1).eq.4.and.iabs(lb(i2)).eq.1))then
              if(iabs(lb(i1)).eq.1)then
        ii = i1
       IF(X2.LE.0.33)THEN
       lb(i1)=7
       e(i1)=dm
       lb(i2)=27
       e(i2)=arho
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i1)=8
       e(i1)=dm
       lb(i2)=26
       e(i2)=arho
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i1)=9
       e(i1)=dm
       lb(i2)=25
       e(i2)=arho
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.33)THEN
       lb(i2)=7
       e(i2)=dm
       lb(i1)=27
       e(i1)=arho
       go to 40
       ENDIF
       if(X2.gt.0.33.and.X2.le.0.67)then
       lb(i2)=8
       e(i2)=dm
       lb(i1)=26
       e(i1)=arho
       go to 40
       endif
       if(X2.gt.0.67)then
       lb(i2)=9
       e(i2)=dm
       lb(i1)=25
       e(i1)=arho
       go to 40
       endif
              endif
       endif 
       if( ((lb(i1).eq.2.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.2))
     &        .OR. ((lb(i1).eq.-2.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.-2)) )then
              if(iabs(lb(i1)).eq.2)then
        ii = i1
       IF(X2.LE.0.5)THEN
       lb(i1)=6
       e(i1)=dm
       lb(i2)=26
       e(i2)=arho
       go to 40
       ELSE
       lb(i1)=7
       e(i1)=dm
       lb(i2)=25
       e(i2)=arho
       go to 40
       endif
              else
        ii = i2
       IF(X2.LE.0.5)THEN
       lb(i2)=6
       e(i2)=dm
       lb(i1)=26
       e(i1)=arho
       go to 40
       ELSE
       lb(i2)=7
       e(i2)=dm
       lb(i1)=25
       e(i1)=arho
       go to 40
       endif
              endif
       ENDIF
       if((iabs(lb(i1)).eq.2.and.lb(i2).eq.4).
     &  or.(lb(i1).eq.4.and.iabs(lb(i2)).eq.2))then
              if(iabs(lb(i1)).eq.2)then
        ii = i1
       IF(X2.LE.0.33)THEN
       lb(i1)=7
       e(i1)=dm
       lb(i2)=26
       e(i2)=arho
       go to 40
       endif
       if(x2.gt.0.33.and.x2.le.0.67)then       
       lb(i1)=6
       e(i1)=dm
       lb(i2)=27
       e(i2)=arho
       go to 40
       endif
       if(x2.gt.0.67)then
       lb(i1)=8
       e(i1)=dm
       lb(i2)=25
       e(i2)=arho
       endif
              else
        ii = i2
       IF(X2.LE.0.33)THEN
       lb(i2)=7
       e(i2)=dm
       lb(i1)=26
       e(i1)=arho
       go to 40
       endif
       if(x2.le.0.67.and.x2.gt.0.33)then       
       lb(i2)=6
       e(i2)=dm
       lb(i1)=27
       e(i1)=arho
       go to 40
       endif
       if(x2.gt.0.67)then
       lb(i2)=8
       e(i2)=dm
       lb(i1)=25
       e(i1)=arho
       endif
              endif
       endif
                     Endif
       if(iblock.eq.79)then
       aomega=0.782
       dmax=srt-0.782-0.02
       dm=rmass(dmax,iseed)
       if( ((lb(i1).eq.1.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.1))
     &  .OR.((lb(i1).eq.-1.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.-1)) )then
              if(iabs(lb(i1)).eq.1)then
        ii = i1
       lb(i1)=9
       e(i1)=dm
       lb(i2)=28
       e(i2)=aomega
       go to 40
              else
        ii = i2
       lb(i2)=9
       e(i2)=dm
       lb(i1)=28
       e(i1)=aomega
       go to 40
              endif
       endif
       if( ((lb(i1).eq.1.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.1))
     &        .OR. ((lb(i1).eq.-1.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.-1)) )then
              if(iabs(lb(i1)).eq.1)then
        ii = i1
       lb(i1)=7
       e(i1)=dm
       lb(i2)=28
       e(i2)=aomega
       go to 40
              else
        ii = i2
       lb(i2)=7
       e(i2)=dm
       lb(i1)=28
       e(i1)=aomega
       go to 40
              endif
       endif
       if( ((lb(i1).eq.2.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.2))
     &       .OR. ((lb(i1).eq.-2.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.-2)) )then
              if(iabs(lb(i1)).eq.2)then
        ii = i1
       lb(i1)=8
       e(i1)=dm
       lb(i2)=28
       e(i2)=aomega
       go to 40
              else
        ii = i2
       lb(i2)=8
       e(i2)=dm
       lb(i1)=28
       e(i1)=aomega
       go to 40
              endif
       endif
       if((iabs(lb(i1)).eq.1.and.lb(i2).eq.4).
     &  or.(lb(i1).eq.4.and.iabs(lb(i2)).eq.1))then
              if(iabs(lb(i1)).eq.1)then
        ii = i1
       lb(i1)=8
       e(i1)=dm
       lb(i2)=28
       e(i2)=aomega
       go to 40
              else
        ii = i2
       lb(i2)=8
       e(i2)=dm
       lb(i1)=28
       e(i1)=aomega
       go to 40
              endif
       endif 
       if( ((lb(i1).eq.2.and.lb(i2).eq.3).
     &  or.(lb(i1).eq.3.and.lb(i2).eq.2))
     &        .OR. ((lb(i1).eq.-2.and.lb(i2).eq.5).
     &  or.(lb(i1).eq.5.and.lb(i2).eq.-2)) )then
              if(iabs(lb(i1)).eq.2)then
        ii = i1
       lb(i1)=6
       e(i1)=dm
       lb(i2)=28
       e(i2)=aomega
       go to 40
              ELSE
        ii = i2
       lb(i2)=6
       e(i2)=dm
       lb(i1)=28
       e(i1)=aomega
              endif
       ENDIF
       if((iabs(lb(i1)).eq.2.and.lb(i2).eq.4).
     &  or.(lb(i1).eq.4.and.iabs(lb(i2)).eq.2))then
              if(iabs(lb(i1)).eq.2)then
        ii = i1
       lb(i1)=7
       e(i1)=dm
       lb(i2)=28
       e(i2)=aomega
       go to 40
              else
        ii = i2
       lb(i2)=7
       e(i2)=dm
       lb(i1)=26
       e(i1)=arho
       go to 40
              endif
       endif
                     Endif
40       em1=e(i1)
       em2=e(i2)
       if(ianti.eq.1 .and. lb(i1).ge.1 .and. lb(i2).ge.1)then
         lb(ii) = -lb(ii)
           jj = i2
          if(ii .eq. i2)jj = i1
         if(iblock .eq. 77)then
          if(lb(jj).eq.3)then
           lb(jj) = 5
          elseif(lb(jj).eq.5)then
           lb(jj) = 3
          endif
         elseif(iblock .eq. 78)then
          if(lb(jj).eq.25)then
           lb(jj) = 27
          elseif(lb(jj).eq.27)then
           lb(jj) = 25
          endif
         endif
       endif
           endif
50          PR2   = (SRT**2 - EM1**2 - EM2**2)**2
     1                - 4.0 * (EM1*EM2)**2
          IF(PR2.LE.0.)PR2=0.00000001
          PR=SQRT(PR2)/(2.*SRT)
          xptr=0.33*pr
         cc1=ptr(xptr,iseed)
         scheck=pr**2-cc1**2
         if(scheck.lt.0) then
            write(99,*) 'scheck36: ', scheck
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
