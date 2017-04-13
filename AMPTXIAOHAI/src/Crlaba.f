        SUBROUTINE Crlaba(PX,PY,PZ,SRT,brel,brsgm,
     &                        I1,I2,nt,IBLOCK,nchrg,icase)
        PARAMETER (MAXSTR=150001, MAXR=1, AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,
     2  AP2=0.13957,AM0=1.232,PI=3.1415926,CUTOFF=1.8966,AVMASS=0.9383)
        PARAMETER  (AKA=0.498,ALA=1.1157,ASA=1.1974)
        PARAMETER  (ETAM=0.5475, AOMEGA=0.782, ARHO=0.77)
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
      if(icase .eq. 3)then
         rrr=RANART(NSEED)
         if(rrr.lt.brel) then
            IBLOCK=8
        else 
            IBLOCK=100
            if(rrr.lt.(brel+brsgm)) then
               LB(i1) = -15 - int(3 * RANART(NSEED))
               e(i1)=asa
            else
               LB(i1)= -14  
               e(i1)=ala
            endif
            LB(i2) = 3 + int(3 * RANART(NSEED))
            e(i2)=0.138
        endif
      endif
      if(icase .eq. 4)then
         rrr=RANART(NSEED)
         if(rrr.lt.brel) then
            IBLOCK=8
         else    
            IBLOCK=102
            LB(i1) = 23
            LB(i2) = -1 - int(2 * RANART(NSEED))
            if(nchrg.eq.-2) LB(i2) = -6
            if(nchrg.eq. 1) LB(i2) = -9
            e(i1) = aka
            e(i2) = 0.938
            if(nchrg.eq.-2.or.nchrg.eq.1) e(i2)=1.232
         endif
      endif
      EM1=E(I1)
      EM2=E(I2)
      PR2   = (SRT**2 - EM1**2 - EM2**2)**2
     1     - 4.0 * (EM1*EM2)**2
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
      CALL ROTATE(PX0,PY0,PZ0,PX,PY,PZ) 
      RETURN
      END
