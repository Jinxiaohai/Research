      SUBROUTINE crkkpi(I1,I2,XSK1, XSK2, XSK3, XSK4,
     &             XSK5, XSK6, XSK7, XSK8, XSK9, XSK10, XSK11, SIGK,
     &             IBLOCK,lbp1,lbp2,emm1,emm2)
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,AMRHO=0.769,AMOMGA=0.782,
     &  AMETA = 0.5473,
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
       IBLOCK=1907
        X1 = RANART(NSEED) * SIGK
        XSK2 = XSK1 + XSK2
        XSK3 = XSK2 + XSK3
        XSK4 = XSK3 + XSK4
        XSK5 = XSK4 + XSK5
        XSK6 = XSK5 + XSK6
        XSK7 = XSK6 + XSK7
        XSK8 = XSK7 + XSK8
        XSK9 = XSK8 + XSK9
        XSK10 = XSK9 + XSK10
        IF (X1 .LE. XSK1) THEN
           LB(I1) = 3 + int(3 * RANART(NSEED))
           LB(I2) = 3 + int(3 * RANART(NSEED))
           E(I1) = AP2
           E(I2) = AP2
           GOTO 100
        ELSE IF (X1 .LE. XSK2) THEN
           LB(I1) = 3 + int(3 * RANART(NSEED))
           LB(I2) = 25 + int(3 * RANART(NSEED))
           E(I1) = AP2
           E(I2) = AMRHO
           GOTO 100
        ELSE IF (X1 .LE. XSK3) THEN
           LB(I1) = 3 + int(3 * RANART(NSEED))
           LB(I2) = 28
           E(I1) = AP2
           E(I2) = AMOMGA
           GOTO 100
        ELSE IF (X1 .LE. XSK4) THEN
           LB(I1) = 3 + int(3 * RANART(NSEED))
           LB(I2) = 0
           E(I1) = AP2
           E(I2) = AMETA
           GOTO 100
        ELSE IF (X1 .LE. XSK5) THEN
           LB(I1) = 25 + int(3 * RANART(NSEED))
           LB(I2) = 25 + int(3 * RANART(NSEED))
           E(I1) = AMRHO
           E(I2) = AMRHO
           GOTO 100
        ELSE IF (X1 .LE. XSK6) THEN
           LB(I1) = 25 + int(3 * RANART(NSEED))
           LB(I2) = 28
           E(I1) = AMRHO
           E(I2) = AMOMGA
           GOTO 100
        ELSE IF (X1 .LE. XSK7) THEN
           LB(I1) = 25 + int(3 * RANART(NSEED))
           LB(I2) = 0
           E(I1) = AMRHO
           E(I2) = AMETA
           GOTO 100
        ELSE IF (X1 .LE. XSK8) THEN
           LB(I1) = 28
           LB(I2) = 28
           E(I1) = AMOMGA
           E(I2) = AMOMGA
           GOTO 100
        ELSE IF (X1 .LE. XSK9) THEN
           LB(I1) = 28
           LB(I2) = 0
           E(I1) = AMOMGA
           E(I2) = AMETA
           GOTO 100
        ELSE IF (X1 .LE. XSK10) THEN
           LB(I1) = 0
           LB(I2) = 0
           E(I1) = AMETA
           E(I2) = AMETA
        ELSE
          iblock = 222
          call rhores(i1,i2)
          lb(i1) = 29
          e(i2)=0.
        END IF
 100    CONTINUE
        lbp1=lb(i1)
        lbp2=lb(i2)
        emm1=e(i1)
        emm2=e(i2)
      RETURN
      END
