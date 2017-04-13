      SUBROUTINE CRLAN(PX,PY,PZ,SRT,I1,I2,IBLOCK)
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,
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
        IBLOCK=71
        NTAG=0
       if( (lb(i1).ge.14.and.lb(i1).le.17) .OR.
     &     (lb(i2).ge.14.and.lb(i2).le.17) )then
        LB(I1)=21
       else
        LB(I1)=23
       endif
        LB(I2)= 3 + int(3 * RANART(NSEED))
        E(I1)=AKA
        E(I2)=0.138
        EM1=E(I1)
        EM2=E(I2)
        PR2   = (SRT**2 - EM1**2 - EM2**2)**2
     1                - 4.0 * (EM1*EM2)**2
          IF(PR2.LE.0.)PR2=1.e-09
          PR=SQRT(PR2)/(2.*SRT)
          C1   = 1.0 - 2.0 * RANART(NSEED)
          T1   = 2.0 * PI * RANART(NSEED)
      S1   = SQRT( 1.0 - C1**2 )
      CT1  = COS(T1)
      ST1  = SIN(T1)
      PZ   = PR * C1
      PX   = PR * S1*CT1
      PY   = PR * S1*ST1
      RETURN
      END
