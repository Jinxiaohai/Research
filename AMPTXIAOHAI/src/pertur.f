      SUBROUTINE pertur(PX,PY,PZ,SRT,IRUN,I1,I2,nt,kp,icont)
      PARAMETER      (MAXSTR=150001,MAXR=1,PI=3.1415926)
      parameter      (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
      PARAMETER (AMN=0.939457,AMP=0.93828,AP1=0.13496,AP2=0.13957)
      PARAMETER      (AKA=0.498,ALA=1.1157,ASA=1.1974,aks=0.895)
      PARAMETER      (ACAS=1.3213,AOME=1.6724,AMRHO=0.769,AMOMGA=0.782)
      PARAMETER      (AETA=0.548,ADIOMG=3.2288)
      parameter            (maxx=20,maxz=24)
      COMMON   /AA/  R(3,MAXSTR)
      COMMON   /BB/  P(3,MAXSTR)
      COMMON   /CC/  E(MAXSTR)
      COMMON   /EE/  ID(MAXSTR),LB(MAXSTR)
      COMMON   /HH/  PROPER(MAXSTR)
      common /ff/f(-mx:mx,-my:my,-mz:mz,-mpx:mpx,-mpy:mpy,-mpz:mpzp)
      common   /gg/  dx,dy,dz,dpx,dpy,dpz
      COMMON   /INPUT/ NSTAR,NDIRCT,DIR
      COMMON   /NN/NNN
      COMMON   /PA/RPION(3,MAXSTR,MAXR)
      COMMON   /PB/PPION(3,MAXSTR,MAXR)
      COMMON   /PC/EPION(MAXSTR,MAXR)
      COMMON   /PD/LPION(MAXSTR,MAXR)
      COMMON   /PE/PROPI(MAXSTR,MAXR)
      COMMON   /RR/  MASSR(0:MAXR)
      COMMON   /BG/BETAX,BETAY,BETAZ,GAMMA
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      COMMON/RNDF77/NSEED
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
      SAVE   
      px0 = px
      py0 = py
      pz0 = pz
      LB1 = LB(I1)
      EM1 = E(I1)
      X1  = R(1,I1)
      Y1  = R(2,I1)
      Z1  = R(3,I1)
      prob1 = PROPER(I1)
      LB2 = LB(I2)
      EM2 = E(I2)
      X2  = R(1,I2)
      Y2  = R(2,I2)
      Z2  = R(3,I2)
      prob2 = PROPER(I2)
      icont = 1
      icsbel = -1
       if( (lb1.eq.21.or.lb1.eq.23.or.iabs(lb1).eq.30) .and.
     &     (iabs(lb2).ge.14.and.iabs(lb2).le.17) )go to 60
       if( (lb2.eq.21.or.lb2.eq.23.or.iabs(lb2).eq.30) .and.
     &     (iabs(lb1).ge.14.and.iabs(lb1).le.17) )go to 60
        if( (lb1.eq.21.or.lb1.eq.23.or.iabs(lb1).eq.30) .and.
     &      (iabs(lb2).eq.40.or.iabs(lb2).eq.41) )go to 70
        if( (lb2.eq.21.or.lb2.eq.23.or.iabs(lb2).eq.30) .and.
     &      (iabs(lb1).eq.40.or.iabs(lb1).eq.41) )go to 70
       if( (((lb1.ge.3.and.lb1.le.5).or.lb1.eq.0) 
     &        .and.(iabs(lb2).eq.40.or.iabs(lb2).eq.41))
     & .OR. (((lb2.ge.3.and.lb2.le.5).or.lb2.eq.0) 
     &        .and.(iabs(lb1).eq.40.or.iabs(lb1).eq.41)) )go to 90
       if( ((lb1.ge.3.and.lb1.le.5).and.iabs(lb2).eq.45)
     & .OR.((lb2.ge.3.and.lb2.le.5).and.iabs(lb1).eq.45) )go to 110
60         if(iabs(lb1).ge.14 .and. iabs(lb1).le.17)then 
             asap = e(i1)
             akap = e(i2)
             idp = i1
           else
             asap = e(i2)
             akap = e(i1)
             idp = i2
           endif
          app = 0.138
         if(srt .lt. (acas+app))return
          srrt = srt - (acas+app) + (amn+akap)
          pkaon = sqrt(((srrt**2-(amn**2+akap**2))/2./amn)**2 - akap**2)
          sigca = 1.5*( akNPsg(pkaon)+akNPsg(pkaon) )
          pii = sqrt((srt**2-(amn+akap)**2)*(srt**2-(amn-akap)**2))
          pff = sqrt((srt**2-(asap+app)**2)*(srt**2-(asap-app)**2))
         cmat = sigca*pii/pff
         sigpi = cmat*
     &            sqrt((srt**2-(acas+app)**2)*(srt**2-(acas-app)**2))/
     &            sqrt((srt**2-(asap+akap)**2)*(srt**2-(asap-akap)**2))
         sigeta = 0.
        if(srt .gt. (acas+aeta))then
           srrt = srt - (acas+aeta) + (amn+akap)
         pkaon = sqrt(((srrt**2-(amn**2+akap**2))/2./amn)**2 - akap**2)
            sigca = 1.5*( akNPsg(pkaon)+akNPsg(pkaon) )
         cmat = sigca*pii/pff
         sigeta = cmat*
     &            sqrt((srt**2-(acas+aeta)**2)*(srt**2-(acas-aeta)**2))/
     &            sqrt((srt**2-(asap+akap)**2)*(srt**2-(asap-akap)**2))
        endif
         sigca = sigpi + sigeta
         sigpe = 0.
           sig = amax1(sigpe,sigca)     
         ds = sqrt(sig/31.4)
         dsr = ds + 0.1
         ec = (em1+em2+0.02)**2
         call distce(i1,i2,dsr,ds,dt,ec,srt,ic,px,py,pz)
           if(ic .eq. -1)return
          brpp = sigca/sig
          if( (lb1.ge.14.and.lb1.le.17) .or.
     &          (lb2.ge.14.and.lb2.le.17) )then
            lbpp1 = 40 + int(2*RANART(NSEED))
          else
            lbpp1 = -40 - int(2*RANART(NSEED))
          endif
              empp1 = acas
           if(RANART(NSEED) .lt. sigpi/sigca)then
            lbpp2 = 3 + int(3*RANART(NSEED))
            empp2 = 0.138
           else
            lbpp2 = 0
            empp2 = aeta
           endif        
          if(RANART(NSEED) .lt. brpp)then
            icont = 0
            lb(i1) = lbpp1
            e(i1) = empp1
            proper(i1) = brpp
            lb(i2) = lbpp2
            e(i2) = empp2
            proper(i2) = 1.
           endif
             go to 700
70         if(iabs(lb1).eq.40 .or. iabs(lb1).eq.41)then 
             acap = e(i1)
             akap = e(i2)
             idp = i1
           else
             acap = e(i2)
             akap = e(i1)
             idp = i2
           endif
           app = 0.138
           ames = 0.138
         if(srt .lt. (aome+ames))return 
          srrt = srt - (aome+ames) + (amn+akap)
         pkaon = sqrt(((srrt**2-(amn**2+akap**2))/2./amn)**2 - akap**2)
         sigomm = 1.5*( akNPsg(pkaon)+akNPsg(pkaon) )
         cmat = sigomm*
     &          sqrt((srt**2-(amn+akap)**2)*(srt**2-(amn-akap)**2))/
     &          sqrt((srt**2-(asa+app)**2)*(srt**2-(asa-app)**2))
        sigom = cmat*
     &           sqrt((srt**2-(aome+ames)**2)*(srt**2-(aome-ames)**2))/
     &           sqrt((srt**2-(acap+akap)**2)*(srt**2-(acap-akap)**2))
          sigpe = 0.
          sig = amax1(sigpe,sigom)     
         ds = sqrt(sig/31.4)
         dsr = ds + 0.1
         ec = (em1+em2+0.02)**2
         call distce(i1,i2,dsr,ds,dt,ec,srt,ic,px,py,pz)
           if(ic .eq. -1)return
           brpp = sigom/sig
           if( (lb1.ge.40.and.lb1.le.41) .or.
     &           (lb2.ge.40.and.lb2.le.41) )then
            lbpp1 = 45
           else
            lbpp1 = -45
           endif
           empp1 = aome
           lbpp2 = 3 + int(3*RANART(NSEED))
           empp2 = ames
           xrand=RANART(NSEED)
         if(xrand .lt. (proper(idp)*brpp))then
            icont = 0
            lb(i1) = lbpp1
            e(i1) = empp1
            proper(i1) = proper(idp)*brpp
            lb(i2) = lbpp2
            e(i2) = empp2
            proper(i2) = 1.
          elseif(xrand.lt.brpp) then
             e(idp) = 0.
          endif
             go to 700
90         if(iabs(lb1).eq.40 .or. iabs(lb1).eq.41)then 
             acap = e(i1)
             app = e(i2)
             idp = i1
             idn = i2
           else
             acap = e(i2)
             app = e(i1)
             idp = i2
             idn = i1
           endif
            akal = aka
         alas = ala
       if(srt .le. (alas+aka))return
           srrt = srt - (acap+app) + (amn+aka)
         pkaon = sqrt(((srrt**2-(amn**2+aka**2))/2./amn)**2 - aka**2)
         sigca = 1.5*( akNPsg(pkaon)+akNPsg(pkaon) )
         cmat = sigca*
     &          sqrt((srt**2-(amn+aka)**2)*(srt**2-(amn-aka)**2))/
     &          sqrt((srt**2-(alas+0.138)**2)*(srt**2-(alas-0.138)**2))
         sigca = cmat*
     &            sqrt((srt**2-(acap+app)**2)*(srt**2-(acap-app)**2))/
     &            sqrt((srt**2-(alas+aka)**2)*(srt**2-(alas-aka)**2))
            dfr = 1./3.
           if(lb(idn).eq.0)dfr = 1.
        sigcal = sigca*dfr*(srt**2-(alas+aka)**2)*
     &           (srt**2-(alas-aka)**2)/(srt**2-(acap+app)**2)/
     &           (srt**2-(acap-app)**2)
          alas = ASA
       if(srt .le. (alas+aka))then
         sigcas = 0.
       else
           srrt = srt - (acap+app) + (amn+aka)
        pkaon = sqrt(((srrt**2-(amn**2+aka**2))/2./amn)**2 - aka**2)
          sigca = 1.5*( akNPsg(pkaon)+akNPsg(pkaon) )
         cmat = sigca*
     &          sqrt((srt**2-(amn+aka)**2)*(srt**2-(amn-aka)**2))/
     &          sqrt((srt**2-(alas+0.138)**2)*(srt**2-(alas-0.138)**2))
         sigca = cmat*
     &            sqrt((srt**2-(acap+app)**2)*(srt**2-(acap-app)**2))/
     &            sqrt((srt**2-(alas+aka)**2)*(srt**2-(alas-aka)**2))
            dfr = 1.
           if(lb(idn).eq.0)dfr = 3.
        sigcas = sigca*dfr*(srt**2-(alas+aka)**2)*
     &           (srt**2-(alas-aka)**2)/(srt**2-(acap+app)**2)/
     &           (srt**2-(acap-app)**2)
       endif
         sig = sigcal + sigcas
         brpp = 1.                                                   
         ds = sqrt(sig/31.4)
         dsr = ds + 0.1
         ec = (em1+em2+0.02)**2
         call distce(i1,i2,dsr,ds,dt,ec,srt,ic,px,py,pz)
       if(ic .eq. -1)then
         ds = sqrt(20.0/31.4)
         dsr = ds + 0.1
         call distce(i1,i2,dsr,ds,dt,ec,srt,icsbel,px,py,pz)
           if(icsbel .eq. -1)return
            empp1 = EM1
            empp2 = EM2
             go to 700
       endif
       IF(sigcal/sig .GT. RANART(NSEED))THEN  
          if(lb1.eq.40.or.lb1.eq.41.or.lb2.eq.40.or.lb2.eq.41)then
          lbpp1 = 21
           lbpp2 = 14
          else
           lbpp1 = 23
           lbpp2 = -14
          endif
         alas = ala
       ELSE
          if(lb1.eq.40.or.lb1.eq.41.or.lb2.eq.40.or.lb2.eq.41)then
           lbpp1 = 21
            lbpp2 = 15 + int(3 * RANART(NSEED))
          else
            lbpp1 = 23
            lbpp2 = -15 - int(3 * RANART(NSEED))
          endif
         alas = ASA       
        ENDIF
             empp1 = aka  
             empp2 = alas 
          if(RANART(NSEED) .lt. proper(idp))then
            icont = 0
            lb(i1) = lbpp1
            e(i1) = empp1
            proper(i1) = 1.
            lb(i2) = lbpp2
            e(i2) = empp2
            proper(i2) = 1.
             go to 700
           else
            e(idp) = 0.
           endif
          return
110         if(lb1 .eq. 45 .or. lb1 .eq. -45)then 
             aomp = e(i1)
             app = e(i2)
             idp = i1
             idn = i2
           else
             aomp = e(i2)
             app = e(i1)
             idp = i2
             idn = i1
           endif
            akal = aka
       if(srt .le. (acas+aka))return
           srrt = srt - (aome+app) + (amn+aka)
         pkaon = sqrt(((srrt**2-(amn**2+aka**2))/2./amn)**2 - aka**2)
           sigca = 1.5*( akNPsg(pkaon)+akNPsg(pkaon) )
         cmat = sigca*
     &          sqrt((srt**2-(amn+aka)**2)*(srt**2-(amn-aka)**2))/
     &          sqrt((srt**2-(asa+0.138)**2)*(srt**2-(asa-0.138)**2))
         sigom = cmat*
     &            sqrt((srt**2-(aomp+app)**2)*(srt**2-(aomp-app)**2))/
     &            sqrt((srt**2-(acas+aka)**2)*(srt**2-(acas-aka)**2))
           dfr = 2./3.
        sigom = sigom*dfr*(srt**2-(acas+aka)**2)*
     &           (srt**2-(acas-aka)**2)/(srt**2-(aomp+app)**2)/
     &           (srt**2-(aomp-app)**2)
         brpp = 1.
         ds = sqrt(sigom/31.4)
         dsr = ds + 0.1
         ec = (em1+em2+0.02)**2
         call distce(i1,i2,dsr,ds,dt,ec,srt,ic,px,py,pz)
       if(ic .eq. -1)then
         ds = sqrt(20.0/31.4)
         dsr = ds + 0.1
         call distce(i1,i2,dsr,ds,dt,ec,srt,icsbel,px,py,pz)
           if(icsbel .eq. -1)return
            empp1 = EM1
            empp2 = EM2
             go to 700
       endif
           if(lb1.eq.45 .or. lb2.eq.45)then
             lbpp1 = 40 + int(2*RANART(NSEED))
             lbpp2 = 21
            else
            lbpp1 = -40 - int(2*RANART(NSEED))
            lbpp2 = 23
           endif
             empp1 = acas
             empp2 = aka  
          if(RANART(NSEED) .lt. proper(idp))then
            icont = 0
            lb(i1) = lbpp1
            e(i1) = empp1
            proper(i1) = proper(idp)
            lb(i2) = lbpp2
            e(i2) = empp2
            proper(i2) = 1.
           else
            e(idp) = 0.
           endif
             go to 700
700    continue
          PR2   = (SRT**2 - EMpp1**2 - EMpp2**2)**2
     &                - 4.0 * (EMpp1*EMpp2)**2
          IF(PR2.LE.0.)PR2=0.00000001
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
       if(icont .eq. 0)return
              E1CM    = SQRT (EMpp1**2 + PX**2 + PY**2 + PZ**2)
              P1BETA  = PX*BETAX + PY*BETAY + PZ*BETAZ
              TRANSF  = GAMMA * ( GAMMA * P1BETA / (GAMMA + 1) + E1CM )
              Ppt11 = BETAX * TRANSF + PX
              Ppt12 = BETAY * TRANSF + PY
              Ppt13 = BETAZ * TRANSF + PZ
         if(icsbel .ne. -1)then
              p(1,i1) = Ppt11
              p(2,i1) = Ppt12
              p(3,i1) = Ppt13
              E2CM    = SQRT (EMpp2**2 + PX**2 + PY**2 + PZ**2)
              TRANSF  = GAMMA * ( -GAMMA * P1BETA / (GAMMA + 1) + E2CM )
              Ppt21 = BETAX * TRANSF - PX
              Ppt22 = BETAY * TRANSF - PY
              Ppt23 = BETAZ * TRANSF - PZ
              p(1,i2) = Ppt21
              p(2,i2) = Ppt22
              p(3,i2) = Ppt23
             return
          endif
                Xpt=X1
                Ypt=Y1
                Zpt=Z1
               NNN=NNN+1
               PROPI(NNN,IRUN)= proper(idp)*brpp
               LPION(NNN,IRUN)= lbpp1
               EPION(NNN,IRUN)= empp1
                RPION(1,NNN,IRUN)=Xpt
                RPION(2,NNN,IRUN)=Ypt
                RPION(3,NNN,IRUN)=Zpt
               PPION(1,NNN,IRUN)=Ppt11
               PPION(2,NNN,IRUN)=Ppt12
               PPION(3,NNN,IRUN)=Ppt13
               dppion(nnn,irun)=dpertp(i1)*dpertp(i2)
            RETURN
            END
