        SUBROUTINE Crksph(PX,PY,PZ,EC,SRT,
     &     emm1,emm2,lbp1,lbp2,I1,I2,ikkg,ikkl,iblock,
     &     icase,srhoks)
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
        sigela=10.
        sigkm=0.
        if((lb1.ge.25.and.lb1.le.28).or.(lb2.ge.25.and.lb2.le.28)) then
           if(iabs(lb1).eq.30.or.iabs(lb2).eq.30) then
              sigkm=srhoks
           elseif((lb1.eq.23.or.lb1.eq.21.or.lb2.eq.23.or.lb2.eq.21)
     1             .and.srt.gt.(ap2+aks)) then
              sigkm=srhoks
           endif
        endif
        if(srt .lt. (aphi+aka)) then
           sig11=0.
           sig22=0.
        else
         if( (iabs(lb1).eq.30.and.(lb2.ge.3.and.lb2.le.5)) .or.
     &       (iabs(lb2).eq.30.and.(lb1.ge.3.and.lb1.le.5)) )then
              dnr =  18.
              ikkl = 0
              IBLOCK = 225
               sig1 = 2047.042  
               sig2 = 1496.692
       elseif((lb1.eq.23.or.lb1.eq.21.and.(lb2.ge.25.and.lb2.le.27)).or.
     &      (lb2.eq.23.or.lb2.eq.21.and.(lb1.ge.25.and.lb1.le.27)) )then
              dnr =  18.
              ikkl = 1
              IBLOCK = 224
               sig1 = 526.702
               sig2 = 1313.960
         elseif( (iabs(lb1).eq.30.and.(lb2.ge.25.and.lb2.le.27)) .or.
     &           (iabs(lb2).eq.30.and.(lb1.ge.25.and.lb1.le.27)) )then
              dnr =  54.
              ikkl = 0
              IBLOCK = 225
               sig1 = 1371.257
               sig2 = 6999.840
         elseif( ((lb1.eq.23.or.lb1.eq.21) .and. lb2.eq.28).or.
     &           ((lb2.eq.23.or.lb2.eq.21) .and. lb1.eq.28) )then
              dnr = 6.
              ikkl = 1
              IBLOCK = 224
               sig1 = 355.429
               sig2 = 440.558
          else
              dnr = 18.
              ikkl = 0
              IBLOCK = 225
               sig1 = 482.292
               sig2 = 1698.903
          endif
            sig11 = 0.
            sig22 = 0.
            scheck=(srt**2-(e(i1)+e(i2))**2)*(srt**2-(e(i1)-e(i2))**2)
            if(scheck.le.0) then
               write(99,*) 'scheck42: ', scheck
               stop
            endif
            pii=sqrt(scheck)
            scheck=(srt**2-(aphi+aka)**2)*(srt**2-(aphi-aka)**2)
            if(scheck.lt.0) then
               write(99,*) 'scheck43: ', scheck
               scheck=0.
            endif
        pff = sqrt(scheck)
          sig11 = sig1*pff/pii*6./dnr/32./pi/srt**2
          if(srt .gt. aphi+aks)then
        pff = sqrt((srt**2-(aphi+aks)**2)*(srt**2-(aphi-aks)**2))
          sig22 = sig2*pff/pii*18./dnr/32./pi/srt**2
           endif
        endif
         sigks=sig11+sig22+sigela+sigkm
        DSkn=SQRT(sigks/PI/10.)
        dsknr=dskn+0.1
        CALL DISTCE(I1,I2,dsknr,DSkn,DT,EC,SRT,IC,
     1  PX,PY,PZ)
        IF(IC.EQ.-1)return
        icase = 1
        ranx = RANART(NSEED) 
         if(ranx .le. (sigela/sigks))then 
            lbp1=lb1
            emm1=e(i1)
            lbp2=lb2
            emm2=e(i2)
            iblock=111
         elseif(ranx .le. ((sigela+sigkm)/sigks))then 
            lbp1=3+int(3*RANART(NSEED))
            emm1=0.14
            if(lb1.eq.23.or.lb2.eq.23) then
               lbp2=30
               emm2=aks
            elseif(lb1.eq.21.or.lb2.eq.21) then
               lbp2=-30
               emm2=aks
            elseif(lb1.eq.30.or.lb2.eq.30) then
               lbp2=23
               emm2=aka
            else
               lbp2=21
               emm2=aka
            endif
            iblock=112
         elseif(ranx .le. ((sigela+sigkm+sig11)/sigks))then 
            lbp2 = 23
            emm2 = aka
            ikkg = 1
            if(lb1.eq.21.or.lb2.eq.21.or.lb1.eq.-30.or.lb2.eq.-30)then
               lbp2=21
               iblock=iblock-100
            endif
            lbp1 = 29
            emm1 = aphi
         else
            lbp2 = 30
            emm2 = aks
            ikkg = 0
            IBLOCK=IBLOCK+2
            if(lb1.eq.21.or.lb2.eq.21.or.lb1.eq.-30.or.lb2.eq.-30)then
               lbp2=-30
               iblock=iblock-100
            endif
            lbp1 = 29
            emm1 = aphi
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
