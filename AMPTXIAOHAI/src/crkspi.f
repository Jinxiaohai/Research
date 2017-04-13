      SUBROUTINE crkspi(I1,I2,XSK1, XSK2, XSK3, XSK4, SIGK,
     & IBLOCK,lbp1,lbp2,emm1,emm2)
        PARAMETER (MAXSTR=150001,MAXR=1)
          PARAMETER (AP1=0.13496,AP2=0.13957,RHOM = 0.770,PI=3.1415926)
        PARAMETER (AETA=0.548,AMOMGA=0.782)
        parameter (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
        COMMON /AA/ R(3,MAXSTR)
        COMMON /BB/ P(3,MAXSTR)
        COMMON /CC/ E(MAXSTR)
        COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
        common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      COMMON/RNDF77/NSEED
      SAVE   
       IBLOCK=466
        X1 = RANART(NSEED) * SIGK
        XSK2 = XSK1 + XSK2
        XSK3 = XSK2 + XSK3
        XSK4 = XSK3 + XSK4
        IF (X1 .LE. XSK1) THEN
           LB(I1) = 3 + int(3 * RANART(NSEED))
           LB(I2) = 25 + int(3 * RANART(NSEED))
           E(I1) = AP2
           E(I2) = rhom
        ELSE IF (X1 .LE. XSK2) THEN
           LB(I1) = 3 + int(3 * RANART(NSEED))
           LB(I2) = 28
           E(I1) = AP2
           E(I2) = AMOMGA
        ELSE IF (X1 .LE. XSK3) THEN
           LB(I1) = 0
           LB(I2) = 25 + int(3 * RANART(NSEED))
           E(I1) = AETA
           E(I2) = rhom
        ELSE
           LB(I1) = 0
           LB(I2) = 28
           E(I1) = AETA
           E(I2) = AMOMGA
        ENDIF
        if(lb(i1).eq.4) E(I1) = AP1
        lbp1=lb(i1)
        lbp2=lb(i2)
        emm1=e(i1)
        emm2=e(i2)
      RETURN
      END
