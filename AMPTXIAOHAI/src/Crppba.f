      SUBROUTINE Crppba(PX,PY,PZ,SRT,I1,I2,IBLOCK)
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,AMRHO=0.769,AMOMGA=0.782,
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
       call pbarfs(srt,npion,iseed)
       nchrg1=3+int(3*RANART(NSEED))
       nchrg2=3+int(3*RANART(NSEED))
      pmass1=ap1
       pmass2=ap1
       if(nchrg1.eq.3.or.nchrg1.eq.5)pmass1=ap2
       if(nchrg2.eq.3.or.nchrg2.eq.5)pmass2=ap2
       IF(NPION.EQ.2)THEN 
       IBLOCK=1902
       LB(I1)=nchrg1
       E(I1)=pmass1
       LB(I2)=nchrg2
       E(I2)=pmass2
       GO TO 50
       ENDIF
       IF(NPION.EQ.3)THEN 
       IBLOCK=1903
       LB(I1)=nchrg1
       E(I1)=pmass1
       LB(I2)=22+nchrg2
            E(I2)=AMRHO
       GO TO 50
       ENDIF
        IF(NPION.EQ.4)THEN 
       IBLOCK=1904
       if(RANART(NSEED).ge.0.5)then
       LB(I1)=22+nchrg1
       E(I1)=AMRHO
       LB(I2)=22+nchrg2
            E(I2)=AMRHO
       else
       LB(I1)=nchrg1
       E(I1)=pmass1
       LB(I2)=28
            E(I2)=AMOMGA
       endif
       GO TO 50
       ENDIF
        IF(NPION.EQ.5)THEN 
       IBLOCK=1905
        LB(I1)=22+nchrg1
       E(I1)=AMRHO
       LB(I2)=28
       E(I2)=AMOMGA
       GO TO 50
       ENDIF
         IF(NPION.EQ.6)THEN 
       IBLOCK=1906
        LB(I1)=28
       E(I1)=AMOMGA
       LB(I2)=28
          E(I2)=AMOMGA
       ENDIF
50    EM1=E(I1)
      EM2=E(I2)
          PR2   = (SRT**2 - EM1**2 - EM2**2)**2
     1                - 4.0 * (EM1*EM2)**2
          IF(PR2.LE.0.)PR2=1.E-08
          PR=SQRT(PR2)/(2.*SRT)
          C1   = 1.0 - 2.0 * RANART(NSEED)
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
