      SUBROUTINE Crkhyp(PX,PY,PZ,SRT,I1,I2,
     &     XKY1, XKY2, XKY3, XKY4, XKY5,
     &     XKY6, XKY7, XKY8, XKY9, XKY10, XKY11, XKY12, XKY13,
     &     XKY14, XKY15, XKY16, XKY17, SIGK, IKMP,
     &     IBLOCK)
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,AMRHO=0.769,AMOMGA=0.782,APHI=1.02,
     2  AP2=0.13957,AM0=1.232,PI=3.1415926,CUTOFF=1.8966,AVMASS=0.9383)
          parameter (pimass=0.140, AMETA = 0.5473, aka=0.498,
     &     aml=1.116,ams=1.193, AM1440 = 1.44, AM1535 = 1.535)
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
       IBLOCK=1908
        X1 = RANART(NSEED) * SIGK
        XKY2 = XKY1 + XKY2
        XKY3 = XKY2 + XKY3
        XKY4 = XKY3 + XKY4
        XKY5 = XKY4 + XKY5
        XKY6 = XKY5 + XKY6
        XKY7 = XKY6 + XKY7
        XKY8 = XKY7 + XKY8
        XKY9 = XKY8 + XKY9
        XKY10 = XKY9 + XKY10
        XKY11 = XKY10 + XKY11
        XKY12 = XKY11 + XKY12
        XKY13 = XKY12 + XKY13
        XKY14 = XKY13 + XKY14
        XKY15 = XKY14 + XKY15
        XKY16 = XKY15 + XKY16
        IF (X1 .LE. XKY1) THEN
           LB(I1) = 3 + int(3 * RANART(NSEED))
           LB(I2) = 1 + int(2 * RANART(NSEED))
           E(I1) = PIMASS
           E(I2) = AMP
           GOTO 100
        ELSE IF (X1 .LE. XKY2) THEN
           LB(I1) = 3 + int(3 * RANART(NSEED))
           LB(I2) = 6 + int(4 * RANART(NSEED))
           E(I1) = PIMASS
           E(I2) = AM0
           GOTO 100
        ELSE IF (X1 .LE. XKY3) THEN
           LB(I1) = 3 + int(3 * RANART(NSEED))
           LB(I2) = 10 + int(2 * RANART(NSEED))
           E(I1) = PIMASS
           E(I2) = AM1440
           GOTO 100
        ELSE IF (X1 .LE. XKY4) THEN
           LB(I1) = 3 + int(3 * RANART(NSEED))
           LB(I2) = 12 + int(2 * RANART(NSEED))
           E(I1) = PIMASS
           E(I2) = AM1535
           GOTO 100
        ELSE IF (X1 .LE. XKY5) THEN
           LB(I1) = 25 + int(3 * RANART(NSEED))
           LB(I2) = 1 + int(2 * RANART(NSEED))
           E(I1) = AMRHO
           E(I2) = AMP
           GOTO 100
        ELSE IF (X1 .LE. XKY6) THEN
           LB(I1) = 25 + int(3 * RANART(NSEED))
           LB(I2) = 6 + int(4 * RANART(NSEED))
           E(I1) = AMRHO
           E(I2) = AM0
           GOTO 100
        ELSE IF (X1 .LE. XKY7) THEN
           LB(I1) = 25 + int(3 * RANART(NSEED))
           LB(I2) = 10 + int(2 * RANART(NSEED))
           E(I1) = AMRHO
           E(I2) = AM1440
           GOTO 100
        ELSE IF (X1 .LE. XKY8) THEN
           LB(I1) = 25 + int(3 * RANART(NSEED))
           LB(I2) = 12 + int(2 * RANART(NSEED))
           E(I1) = AMRHO
           E(I2) = AM1535
           GOTO 100
        ELSE IF (X1 .LE. XKY9) THEN
           LB(I1) = 28
           LB(I2) = 1 + int(2 * RANART(NSEED))
           E(I1) = AMOMGA
           E(I2) = AMP
           GOTO 100
        ELSE IF (X1 .LE. XKY10) THEN
           LB(I1) = 28
           LB(I2) = 6 + int(4 * RANART(NSEED))
           E(I1) = AMOMGA
           E(I2) = AM0
           GOTO 100
        ELSE IF (X1 .LE. XKY11) THEN
           LB(I1) = 28
           LB(I2) = 10 + int(2 * RANART(NSEED))
           E(I1) = AMOMGA
           E(I2) = AM1440
           GOTO 100
        ELSE IF (X1 .LE. XKY12) THEN
           LB(I1) = 28
           LB(I2) = 12 + int(2 * RANART(NSEED))
           E(I1) = AMOMGA
           E(I2) = AM1535
           GOTO 100
        ELSE IF (X1 .LE. XKY13) THEN
           LB(I1) = 0
           LB(I2) = 1 + int(2 * RANART(NSEED))
           E(I1) = AMETA
           E(I2) = AMP
           GOTO 100
        ELSE IF (X1 .LE. XKY14) THEN
           LB(I1) = 0
           LB(I2) = 6 + int(4 * RANART(NSEED))
           E(I1) = AMETA
           E(I2) = AM0
           GOTO 100
        ELSE IF (X1 .LE. XKY15) THEN
           LB(I1) = 0
           LB(I2) = 10 + int(2 * RANART(NSEED))
           E(I1) = AMETA
           E(I2) = AM1440
           GOTO 100
        ELSE IF (X1 .LE. XKY16) THEN
           LB(I1) = 0
           LB(I2) = 12 + int(2 * RANART(NSEED))
           E(I1) = AMETA
           E(I2) = AM1535
           GOTO 100
        ELSE
           LB(I1) = 29
           LB(I2) = 1 + int(2 * RANART(NSEED))
           E(I1) = APHI
           E(I2) = AMN
          IBLOCK=222
           GOTO 100
        END IF
 100    CONTINUE
         if(IKMP .eq. -1) LB(I2) = -LB(I2)
      EM1=E(I1)
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
