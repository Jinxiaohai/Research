      SUBROUTINE CREN(PX,PY,PZ,SRT,I1,I2,IBLOCK)
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
        NTAG=0
        IBLOCK=7
        ianti=0
        if(lb(i1).lt.0 .or. lb(i2).lt.0)then
          ianti=1
          iblock=-7
        endif
       KAONC=0
       IF(PNLKA(SRT)/(PNLKA(SRT)
     & +PNSKA(SRT)).GT.RANART(NSEED))KAONC=1
       IF(E(I1).LE.0.6)THEN
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
            lb(i1)=21
            lb(i2)=-lb(i2)
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
            lb(i2)=21
            lb(i1)=-lb(i1)
          endif
       ENDIF
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
