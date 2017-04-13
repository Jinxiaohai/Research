        SUBROUTINE Crkphi(PX,PY,PZ,EC,SRT,IBLOCK,
     &                  emm1,emm2,lbp1,lbp2,I1,I2,ikk,icase,rrkk,prkk)
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,AP2=0.13957,APHI=1.02,
     2  AM0=1.232,AMNS=1.52,PI=3.1415926,CUTOFF=1.8966,AVMASS=0.9383)
        PARAMETER      (AKA=0.498,ALA=1.1157,ASA=1.1974,ACAS=1.3213)
        PARAMETER      (AKS=0.895,AOMEGA=0.7819, ARHO=0.77)
        parameter     (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
        COMMON /AA/ R(3,MAXSTR)
        COMMON /BB/ P(3,MAXSTR)
        COMMON /CC/ E(MAXSTR)
        COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
        common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      COMMON/RNDF77/NSEED
      SAVE   
        lb1 = lb(i1) 
        lb2 = lb(i2) 
        icase = 0
        if(srt .lt. (aphi+ap1)) then
           sig1 = 0.
           sig2 = 0.
           sig3 = 0.
        else
         if((lb1.eq.23.and.lb2.eq.21).or.(lb2.eq.23.and.lb1.eq.21))then
            dnr =  4.
            ikk = 2
          elseif((lb1.eq.21.and.lb2.eq.30).or.(lb2.eq.21.and.lb1.eq.30)
     & .or.(lb1.eq.23.and.lb2.eq.-30).or.(lb2.eq.23.and.lb1.eq.-30))then
             dnr = 12.
             ikk = 1
          else
             dnr = 36.
             ikk = 0
          endif
          sig1 = 0.
          sig2 = 0.
          sig3 = 0.
          srri = E(i1)+E(i2)
          srr1 = aphi+ap1
          srr2 = aphi+aomega
          srr3 = aphi+arho
          pii = (srt**2-(e(i1)+e(i2))**2)*(srt**2-(e(i1)-e(i2))**2)
          srrt = srt - amax1(srri,srr1)
          if(srrt .lt. 0.3 .and. srrt .gt. 0.01)then
          sig = 1.69/(srrt**0.141 - 0.407)
         else
          sig = 3.74 + 0.008*srrt**1.9
         endif                 
          sig1=sig*(9./dnr)*(srt**2-(aphi+ap1)**2)*
     &           (srt**2-(aphi-ap1)**2)/pii
          if(srt .gt. aphi+aomega)then
          srrt = srt - amax1(srri,srr2)
          if(srrt .lt. 0.3 .and. srrt .gt. 0.01)then
          sig = 1.69/(srrt**0.141 - 0.407)
         else
          sig = 3.74 + 0.008*srrt**1.9
         endif                 
          sig2=sig*(9./dnr)*(srt**2-(aphi+aomega)**2)*
     &           (srt**2-(aphi-aomega)**2)/pii
           endif
         if(srt .gt. aphi+arho)then
          srrt = srt - amax1(srri,srr3)
          if(srrt .lt. 0.3 .and. srrt .gt. 0.01)then
          sig = 1.69/(srrt**0.141 - 0.407)
         else
          sig = 3.74 + 0.008*srrt**1.9
         endif                 
          sig3=sig*(27./dnr)*(srt**2-(aphi+arho)**2)*
     &           (srt**2-(aphi-arho)**2)/pii
         endif                 
        endif
        rrkk0=rrkk
        prkk0=prkk
        SIGM=0.
        if((lb1.eq.23.and.lb2.eq.21).or.(lb2.eq.23.and.lb1.eq.21))then
           CALL XKKANN(SRT, XSK1, XSK2, XSK3, XSK4, XSK5,
     &          XSK6, XSK7, XSK8, XSK9, XSK10, XSK11, SIGM, rrkk0)
        elseif((lb1.eq.21.and.lb2.eq.30).or.(lb2.eq.21.and.lb1.eq.30)
     & .or.(lb1.eq.23.and.lb2.eq.-30).or.(lb2.eq.23.and.lb1.eq.-30))then
           CALL XKKSAN(i1,i2,SRT,SIGKS1,SIGKS2,SIGKS3,SIGKS4,SIGM,prkk0)
        else
        endif
        sigm0=sigm
        sigks = sig1 + sig2 + sig3 + SIGM
        DSkn=SQRT(sigks/PI/10.)
        dsknr=dskn+0.1
        CALL DISTCE(I1,I2,dsknr,DSkn,DT,EC,SRT,IC,
     1  PX,PY,PZ)
        IF(IC.EQ.-1)return
        icase = 1
        ranx = RANART(NSEED) 
        lbp1 = 29
        emm1 = aphi
        if(ranx .le. sig1/sigks)then 
           lbp2 = 3 + int(3*RANART(NSEED))
           emm2 = ap1
        elseif(ranx .le. (sig1+sig2)/sigks)then
           lbp2 = 28
           emm2 = aomega
        elseif(ranx .le. (sig1+sig2+sig3)/sigks)then
           lbp2 = 25 + int(3*RANART(NSEED))
           emm2 = arho
        else
           if((lb1.eq.23.and.lb2.eq.21)
     &          .or.(lb2.eq.23.and.lb1.eq.21))then
              CALL crkkpi(I1,I2,XSK1, XSK2, XSK3, XSK4,
     &             XSK5, XSK6, XSK7, XSK8, XSK9, XSK10, XSK11, SIGM0,
     &             IBLOCK,lbp1,lbp2,emm1,emm2)
           elseif((lb1.eq.21.and.lb2.eq.30)
     &             .or.(lb2.eq.21.and.lb1.eq.30)
     &             .or.(lb1.eq.23.and.lb2.eq.-30)
     &             .or.(lb2.eq.23.and.lb1.eq.-30))then
              CALL crkspi(I1,I2,SIGKS1, SIGKS2, SIGKS3, SIGKS4,
     &             SIGM0,IBLOCK,lbp1,lbp2,emm1,emm2)
           else
           endif
        endif
        PX0=PX
        PY0=PY
        PZ0=PZ
           PR2   = (SRT**2 - EMM1**2 - EMM2**2)**2
     1                - 4.0 * (EMM1*EMM2)**2
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
