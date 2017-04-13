      SUBROUTINE Crdir(PX,PY,PZ,SRT,I1,I2,IBLOCK)
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
        IBLOCK=999
        NTAG=0
        EM1=E(I1)
        EM2=E(I2)
        PR2   = (SRT**2 - EM1**2 - EM2**2)**2
     1                - 4.0 * (EM1*EM2)**2
          IF(PR2.LE.0.)PR2=1.e-09
          PR=SQRT(PR2)/(2.*SRT)
          xptr=0.33*pr
         cc1=ptr(xptr,iseed)
         scheck=pr**2-cc1**2
         if(scheck.lt.0) then
            write(99,*) 'scheck37: ', scheck
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
      call rotate(px0,py0,pz0,px,py,pz)
      RETURN
      END
