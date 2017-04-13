        SUBROUTINE Crkpla(PX,PY,PZ,EC,SRT,spika,
     &                  emm1,emm2,lbp1,lbp2,I1,I2,icase,srhoks)
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,AP2=0.13957,AMRHO=0.769,AMOMGA=0.782,
     2  AM0=1.232,PI=3.1415926,CUTOFF=1.8966,AVMASS=0.9383)
        PARAMETER (AKA=0.498,AKS=0.895,ALA=1.1157,ASA=1.1974
     1 ,APHI=1.02)
        PARAMETER (AM1440 = 1.44, AM1535 = 1.535)
        parameter     (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
        COMMON /AA/ R(3,MAXSTR)
        COMMON /BB/ P(3,MAXSTR)
        COMMON /CC/ E(MAXSTR)
        COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
        common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      COMMON/RNDF77/NSEED
      SAVE   
          emm1=0.
          emm2=0.
          lbp1=0
          lbp2=0
           XKP0 = spika
           XKP1 = 0.
           XKP2 = 0.
           XKP3 = 0.
           XKP4 = 0.
           XKP5 = 0.
           XKP6 = 0.
           XKP7 = 0.
           XKP8 = 0.
           XKP9 = 0.
           XKP10 = 0.
           sigm = 15.
        pdd = (srt**2-(aka+ap1)**2)*(srt**2-(aka-ap1)**2)
         if(srt .lt. (ala+amn))go to 70
        XKP1 = sigm*(4./3.)*(srt**2-(ala+amn)**2)*
     &           (srt**2-(ala-amn)**2)/pdd
         if(srt .gt. (ala+am0))then
        XKP2 = sigm*(16./3.)*(srt**2-(ala+am0)**2)*
     &           (srt**2-(ala-am0)**2)/pdd
         endif
         if(srt .gt. (ala+am1440))then
        XKP3 = sigm*(4./3.)*(srt**2-(ala+am1440)**2)*
     &           (srt**2-(ala-am1440)**2)/pdd
         endif
         if(srt .gt. (ala+am1535))then
        XKP4 = sigm*(4./3.)*(srt**2-(ala+am1535)**2)*
     &           (srt**2-(ala-am1535)**2)/pdd
         endif
         if(srt .gt. (asa+amn))then
        XKP5 = sigm*4.*(srt**2-(asa+amn)**2)*
     &           (srt**2-(asa-amn)**2)/pdd
         endif
         if(srt .gt. (asa+am0))then
        XKP6 = sigm*16.*(srt**2-(asa+am0)**2)*
     &           (srt**2-(asa-am0)**2)/pdd
         endif
         if(srt .gt. (asa+am1440))then
        XKP7 = sigm*4.*(srt**2-(asa+am1440)**2)*
     &           (srt**2-(asa-am1440)**2)/pdd
         endif
         if(srt .gt. (asa+am1535))then
        XKP8 = sigm*4.*(srt**2-(asa+am1535)**2)*
     &           (srt**2-(asa-am1535)**2)/pdd
         endif
70     continue
          sig1 = 195.639
          sig2 = 372.378
       if(srt .gt. aphi+aka)then
        pff = sqrt((srt**2-(aphi+aka)**2)*(srt**2-(aphi-aka)**2))
        scheck=pdd
        if(scheck.le.0) then
           write(99,*) 'scheck40: ', scheck
           stop
        endif
         XKP9 = sig1*pff/sqrt(pdd)*1./32./pi/srt**2
        if(srt .gt. aphi+aks)then
        pff = sqrt((srt**2-(aphi+aks)**2)*(srt**2-(aphi-aks)**2))
        scheck=pdd
        if(scheck.le.0) then
           write(99,*) 'scheck41: ', scheck
           stop
        endif
         XKP10 = sig2*pff/sqrt(pdd)*3./32./pi/srt**2
       endif
        endif
        sigpik=0.
        if(srt.gt.(amrho+aks)) then
           sigpik=srhoks*9.
     1          *(srt**2-(0.77-aks)**2)*(srt**2-(0.77+aks)**2)/4
     2          /srt**2/(px**2+py**2+pz**2)
           if(srt.gt.(amomga+aks)) sigpik=sigpik*12./9.
        endif
         sigkp = XKP0 + XKP1 + XKP2 + XKP3 + XKP4
     &         + XKP5 + XKP6 + XKP7 + XKP8 + XKP9 + XKP10 +sigpik
           icase = 0 
         DSkn=SQRT(sigkp/PI/10.)
        dsknr=dskn+0.1
        CALL DISTCE(I1,I2,dsknr,DSkn,DT,EC,SRT,IC,
     1  PX,PY,PZ)
        IF(IC.EQ.-1)return
        randu = RANART(NSEED)*sigkp
        XKP1 = XKP0 + XKP1
        XKP2 = XKP1 + XKP2
        XKP3 = XKP2 + XKP3
        XKP4 = XKP3 + XKP4
        XKP5 = XKP4 + XKP5
        XKP6 = XKP5 + XKP6
        XKP7 = XKP6 + XKP7
        XKP8 = XKP7 + XKP8
        XKP9 = XKP8 + XKP9
        XKP10 = XKP9 + XKP10
         if(randu .le. XKP0)then
           icase = 1
            return
         else
           icase = 2
         if( randu .le. XKP1 )then
             lbp1 = -14
             lbp2 = 1 + int(2*RANART(NSEED))
             emm1 = ala
             emm2 = amn
             go to 60
         elseif( randu .le. XKP2 )then
             lbp1 = -14
             lbp2 = 6 + int(4*RANART(NSEED))
             emm1 = ala
             emm2 = am0
             go to 60
         elseif( randu .le. XKP3 )then
             lbp1 = -14
             lbp2 = 10 + int(2*RANART(NSEED))
             emm1 = ala
             emm2 = am1440
             go to 60
         elseif( randu .le. XKP4 )then
             lbp1 = -14
             lbp2 = 12 + int(2*RANART(NSEED))
             emm1 = ala
             emm2 = am1535
             go to 60
         elseif( randu .le. XKP5 )then
             lbp1 = -15 - int(3*RANART(NSEED))
             lbp2 = 1 + int(2*RANART(NSEED))
             emm1 = asa
             emm2 = amn
             go to 60
         elseif( randu .le. XKP6 )then
             lbp1 = -15 - int(3*RANART(NSEED))
             lbp2 = 6 + int(4*RANART(NSEED))
             emm1 = asa
             emm2 = am0
             go to 60
          elseif( randu .lt. XKP7 )then
             lbp1 = -15 - int(3*RANART(NSEED))
             lbp2 = 10 + int(2*RANART(NSEED))
             emm1 = asa
             emm2 = am1440
             go to 60
          elseif( randu .lt. XKP8 )then
             lbp1 = -15 - int(3*RANART(NSEED))
             lbp2 = 12 + int(2*RANART(NSEED))
             emm1 = asa
             emm2 = am1535
             go to 60
          elseif( randu .lt. XKP9 )then
            icase = 3
             lbp1 = 29
             lbp2 = 23
             emm1 = aphi
             emm2 = aka
           if(lb(i1).eq.21.or.lb(i2).eq.21)then
             lbp2 = 21
             icase = -3
           endif
             go to 60
          elseif( randu .lt. XKP10 )then
            icase = 4
             lbp1 = 29
             lbp2 = 30
             emm1 = aphi
             emm2 = aks
           if(lb(i1).eq.21.or.lb(i2).eq.21)then
             lbp2 = -30
             icase = -4
           endif
           go to 60
          else
            icase=5
            lbp1=25+int(3*RANART(NSEED))
            lbp2=30
            emm1=amrho
            emm2=aks
            if(srt.gt.(amomga+aks).and.RANART(NSEED).lt.0.25) then
               lbp1=28
               emm1=amomga
            endif
            if(lb(i1).eq.21.or.lb(i2).eq.21)then
               lbp2=-30
               icase=-5
            endif
          endif
          endif
60       if( icase.eq.2 .and. (lb(i1).eq.21.or.lb(i2).eq.21) )then
            lbp1 = -lbp1
            lbp2 = -lbp2
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
