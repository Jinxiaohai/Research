       SUBROUTINE ARTMN
      PARAMETER     (MAXSTR=150001,MAXR=1,AMU= 0.9383,
     1               AKA=0.498,etaM=0.5475)
      PARAMETER     (MAXX   =   20,  MAXZ  =    24)
      PARAMETER     (ISUM   =   1001,  IGAM  =    1100)
      parameter     (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
      INTEGER   OUTPAR, zta,zpr
      COMMON  /AA/      R(3,MAXSTR)
      COMMON  /BB/      P(3,MAXSTR)
      COMMON  /CC/      E(MAXSTR)
      COMMON  /DD/      RHO(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHOP(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHON(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ)
      COMMON  /EE/      ID(MAXSTR),LB(MAXSTR)
      COMMON  /HH/  PROPER(MAXSTR)
      common  /ff/f(-mx:mx,-my:my,-mz:mz,-mpx:mpx,-mpy:mpy,-mpz:mpzp)
      common  /gg/      dx,dy,dz,dpx,dpy,dpz
      COMMON  /INPUT/ NSTAR,NDIRCT,DIR
      COMMON  /PP/      PRHO(-20:20,-24:24)
      COMMON  /QQ/      PHRHO(-MAXZ:MAXZ,-24:24)
      COMMON  /RR/      MASSR(0:MAXR)
      common  /ss/      inout(20)
      common  /zz/      zta,zpr
      COMMON  /RUN/     NUM
      COMMON  /KKK/     TKAON(7),EKAON(7,0:2000)
      COMMON  /KAON/    AK(3,50,36),SPECK(50,36,7),MF
      COMMON/TABLE/ xarray(0:1000),earray(0:1000)
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      COMMON  /DDpi/    piRHO(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ)
      common  /tt/  PEL(-maxx:maxx,-maxx:maxx,-maxz:maxz)
     &,rxy(-maxx:maxx,-maxx:maxx,-maxz:maxz)
      DIMENSION TEMP(3,MAXSTR),SKAON(7),SEKAON(7,0:2000)
      COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &   IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
      COMMON /INPUT3/ PLAB, ELAB, ZEROPT, B0, BI, BM, DENCUT, CYCBOX
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)
      COMMON /ARERCP/PRO1(MAXSTR, MAXR)
      COMMON /ARERC1/MULTI1(MAXR)
      COMMON /ARPRC1/ITYP1(MAXSTR, MAXR),
     &     GX1(MAXSTR, MAXR), GY1(MAXSTR, MAXR), GZ1(MAXSTR, MAXR), 
     &     FT1(MAXSTR, MAXR),
     &     PX1(MAXSTR, MAXR), PY1(MAXSTR, MAXR), PZ1(MAXSTR, MAXR),
     &     EE1(MAXSTR, MAXR), XM1(MAXSTR, MAXR)
      DIMENSION NPI(MAXR)
      DIMENSION RT(3, MAXSTR, MAXR), PT(3, MAXSTR, MAXR)
     &     , ET(MAXSTR, MAXR), LT(MAXSTR, MAXR), PROT(MAXSTR, MAXR)
      EXTERNAL IARFLV, INVFLV
      common /lastt/itimeh,bimp 
      common/snn/efrm,npart1,npart2,epsiPz,epsiPt,PZPROJ,PZTARG
      COMMON/hbt/lblast(MAXSTR),xlast(4,MAXSTR),plast(4,MAXSTR),nlast
      common/resdcy/NSAV,iksdcy
      COMMON/RNDF77/NSEED
      COMMON/FTMAX/ftsv(MAXSTR),ftsvt(MAXSTR, MAXR)
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      real zet(-45:45)
      SAVE   
      data zet /
     4     1.,0.,0.,0.,0.,
     3     1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     2     -1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     1     0.,0.,0.,-1.,0.,1.,0.,-1.,0.,-1.,
     s     0.,-2.,-1.,0.,1.,0.,0.,0.,0.,-1.,
     e     0.,
     s     1.,0.,-1.,0.,1.,-1.,0.,1.,2.,0.,
     1     1.,0.,1.,0.,-1.,0.,1.,0.,0.,0.,
     2     -1.,0.,1.,0.,-1.,0.,1.,0.,0.,1.,
     3     0.,0.,0.,0.,0.,0.,0.,0.,0.,-1.,
     4     0.,0.,0.,0.,-1./
      nlast=0
      do 1002 i=1,MAXSTR
         ftsv(i)=0.
         do 1101 irun=1,maxr
            ftsvt(i,irun)=0.
 1101    continue
         lblast(i)=999
         do 1001 j=1,4
            xlast(j,i)=0.
            plast(j,i)=0.
 1001    continue
 1002 continue
       call tablem
       ikaon=1
       nstar=1
       ndirct=0
       dir=0.02
       asy=0.032
       ESBIN=0.04
       MF=36
      RADTA  = 1.124 * FLOAT(MASSTA)**(1./3.)
      RADPR  = 1.124 * FLOAT(MASSPR)**(1./3.)
      ZDIST  = RADTA + RADPR
      BMAX   = RADTA + RADPR
      MASS   = MASSTA + MASSPR
      NTOTAL = NUM * MASS
      IF (NTOTAL .GT. MAXSTR) THEN
        WRITE(12,'(//10X,''**** FATAL ERROR: TOO MANY TEST PART. ****'//
     & ' '')')
        STOP
      END IF
      ETA    = FLOAT(MASSTA) * AMU
      PZTA   = 0.0
      BETATA = 0.0
      GAMMTA = 1.0
      EPR    = FLOAT(MASSPR) * (AMU + 0.001 * ELAB)
      PZPR   = SQRT( EPR**2 - (AMU * FLOAT(MASSPR))**2 )
      BETAPR = PZPR / EPR
      GAMMPR = 1.0 / SQRT( 1.0 - BETAPR**2 )
        BETAC=(PZPR+PZTA)/(EPR+ETA)
        GAMMC=1.0 / SQRT(1.-BETAC**2)
      IF (INSYS .NE. 0) THEN
        S      = (EPR+ETA)**2 - PZPR**2
        xx1=4.*alog(float(massta))
        xx2=4.*alog(float(masspr))
        xx1=exp(xx1)
        xx2=exp(xx2)
        PSQARE = (S**2 + (xx1+ xx2) * AMU**4
     &             - 2.0 * S * AMU**2 * FLOAT(MASSTA**2 + MASSPR**2)
     &             - 2.0 * FLOAT(MASSTA**2 * MASSPR**2) * AMU**4)
     &           / (4.0 * S)
        ETA    = SQRT ( PSQARE + (FLOAT(MASSTA) * AMU)**2 )
        PZTA   = - SQRT(PSQARE)
        BETATA = PZTA / ETA
        GAMMTA = 1.0 / SQRT( 1.0 - BETATA**2 )
        EPR    = SQRT ( PSQARE + (FLOAT(MASSPR) * AMU)**2 )
        PZPR   = SQRT(PSQARE)
        BETAPR = PZPR/ EPR
        GAMMPR = 1.0 / SQRT( 1.0 - BETAPR**2 )
      ELSE
      END IF
      PZTA = PZTA / FLOAT(MASSTA)
      PZPR = PZPR / FLOAT(MASSPR)
      ECMS0=ETA+EPR
       DO 50000 IMANY=1,MANYB
       if (manyb. gt.1) then
111       BX=1.0-2.0*RANART(NSEED)
       BY=1.0-2.0*RANART(NSEED)
       B2=BX*BX+BY*BY
       IF(B2.GT.1.0) GO TO 111       
       B=SQRT(B2)*(BM-BI)+BI
       ELSE
       B=B0
       ENDIF
      call coulin(masspr,massta,NUM)
      CALL INIT(1       ,MASSTA   ,NUM     ,RADTA,
     &          B/2.    ,ZEROPT+ZDIST/2.   ,PZTA,
     &          GAMMTA  ,ISEED    ,MASS    ,IMOMEN)
      CALL INIT(1+MASSTA,MASS     ,NUM     ,RADPR,
     &          -B/2.   ,ZEROPT-ZDIST/2.   ,PZPR,
     &          GAMMPR  ,ISEED    ,MASS    ,IMOMEN)
      OUTPAR = 0
      MASSR(0)=0
      DO 1003 IR =1,NUM
      MASSR(IR)=MASS
 1003 CONTINUE
      CALL DENS(IPOT,MASS,NUM,OUTPAR)
      IF (ICOLL .NE. -1) THEN
        DO 700 I = 1,NTOTAL
          IX = NINT( R(1,I) )
          IY = NINT( R(2,I) )
          IZ = NINT( R(3,I) )
          IF(IX.GE.MAXX.OR.IY.GE.MAXX.OR.IZ.GE.MAXZ
     1         .OR.IX.LE.-MAXX.OR.IY.LE.-MAXX.OR.IZ.LE.-MAXZ) goto 700
          CALL GRADU(IPOT,IX,IY,IZ,GRADX,GRADY,GRADZ)
          P(1,I) = P(1,I) - (0.5 * DT) * GRADX
          P(2,I) = P(2,I) - (0.5 * DT) * GRADY
          P(3,I) = P(3,I) - (0.5 * DT) * GRADZ
  700   CONTINUE
      END IF
        RCNNE  = 0
       RDD  = 0
       RPP  = 0
       rppk = 0
       RPN  = 0
       rpd  = 0
       RKN  = 0
       RNNK = 0
       RDDK = 0
       RNDK = 0
      RCNND  = 0
      RCNDN  = 0
      RCOLL  = 0
      RBLOC  = 0
      RDIRT  = 0
      RDECAY = 0
      RRES   = 0
      DO 1005 KKK=1,5
         SKAON(KKK)  = 0
         DO 1004 IS=1,2000
            SEKAON(KKK,IS)=0
 1004    CONTINUE
 1005 CONTINUE
       pr0=0.
       pr1=0.
       ska0=0.
       ska1=0.
      IF (IAPAR2(1) .NE. 1) THEN
         DO 1016 I = 1, MAXSTR
            DO 1015 J = 1, 3
               R(J, I) = 0.
               P(J, I) = 0.
 1015       CONTINUE
            E(I) = 0.
            LB(I) = 0
            ID(I)=0
           PROPER(I) = 1.
 1016   CONTINUE
         MASS = 0
         NP = 0
         DO 1017 J = 1, NUM
            MASSR(J) = 0
            NPI(J) = 1
 1017    CONTINUE
         DO 1019 I = 1, MAXR
            DO 1018 J = 1, MAXSTR
               RT(1, J, I) = 0.
               RT(2, J, I) = 0.
               RT(3, J, I) = 0.
               PT(1, J, I) = 0.
               PT(2, J, I) = 0.
               PT(3, J, I) = 0.
               ET(J, I) = 0.
               LT(J, I) = 0
               PROT(J, I) = 1.
 1018       CONTINUE
 1019    CONTINUE
      END IF
      DO 10000 NT = 1,NTMAX
      LP1=0
      LP2=0
      LP3=0
      LD1=0
      LD2=0
      LD3=0
      LD4=0
      LN1=0
      LN2=0
      LN5=0
      LE=0
      LKAON=0
      LKAONS=0
        IF (ICOLL .NE. 1) THEN
           numnt=nt
          CALL RELCOL(LCOLL,LBLOC,LCNNE,LDD,LPP,lppk,
     &    LPN,lpd,LRHO,LOMEGA,LKN,LNNK,LDDK,LNDK,LCNND,
     &    LCNDN,LDIRT,LDECAY,LRES,LDOU,LDDRHO,LNNRHO,
     &    LNNOM,numnt,ntmax,sp,akaon,sk)
          RCOLL = RCOLL + FLOAT(LCOLL)/num
          RBLOC = RBLOC + FLOAT(LBLOC)/num
          RCNNE = RCNNE + FLOAT(LCNNE)/num
         RDD   = RDD   + FLOAT(LDD)/num
         RPP   = RPP   + FLOAT(LPP)/NUM
         rppk  =rppk   + float(lppk)/num
         RPN   = RPN   + FLOAT(LPN)/NUM
         rpd   =rpd    + float(lpd)/num
         RKN   = RKN   + FLOAT(LKN)/NUM
         RNNK  =RNNK   + FLOAT(LNNK)/NUM
         RDDK  =RDDK   + FLOAT(LDDK)/NUM
         RNDK  =RNDK   + FLOAT(LNDK)/NUM
          RCNND = RCNND + FLOAT(LCNND)/num
          RCNDN = RCNDN + FLOAT(LCNDN)/num
          RDIRT = RDIRT + FLOAT(LDIRT)/num
          RDECAY= RDECAY+ FLOAT(LDECAY)/num
          RRES  = RRES  + FLOAT(LRES)/num
          ADIRT=LDIRT/DT/num
          ACOLL=(LCOLL-LBLOC)/DT/num
          ACNND=LCNND/DT/num
          ACNDN=LCNDN/DT/num
          ADECAY=LDECAY/DT/num
          ARES=LRES/DT/num
         ADOU=LDOU/DT/NUM
         ADDRHO=LDDRHO/DT/NUM
         ANNRHO=LNNRHO/DT/NUM
         ANNOM=LNNOM/DT/NUM
         ADD=LDD/DT/num
         APP=LPP/DT/num
         appk=lppk/dt/num
          APN=LPN/DT/num
         apd=lpd/dt/num
         arh=lrho/dt/num
         aom=lomega/dt/num
         AKN=LKN/DT/num
         ANNK=LNNK/DT/num
         ADDK=LDDK/DT/num
         ANDK=LNDK/DT/num
        END IF
        CALL DENS(IPOT,MASS,NUM,OUTPAR)
       sumene=0
        ISO=0
        DO 201 MRUN=1,NUM
        ISO=ISO+MASSR(MRUN-1)
        DO 201 I0=1,MASSR(MRUN)
        I =I0+ISO
        ETOTAL = SQRT( E(I)**2 + P(1,I)**2 + P(2,I)**2 +P(3,I)**2 )
       sumene=sumene+etotal
         if(kpoten.ne.0.and.lb(i).eq.23)then
             den=0.
              IX = NINT( R(1,I) )
              IY = NINT( R(2,I) )
              IZ = NINT( R(3,I) )
              IF(IX.LT.MAXX.AND.IY.LT.MAXX.AND.IZ.LT.MAXZ
     1             .AND.IX.GT.-MAXX.AND.IY.GT.-MAXX.AND.IZ.GT.-MAXZ)
     2             den=rho(ix,iy,iz)
          akg = 0.1727
          bkg = 0.333
         rnsg = den
         ecor = - akg*rnsg + (bkg*den)**2
         etotal = sqrt(etotal**2 + ecor)
         endif
         if(kpoten.ne.0.and.lb(i).eq.21)then
             den=0.
              IX = NINT( R(1,I) )
              IY = NINT( R(2,I) )
              IZ = NINT( R(3,I) )
              IF(IX.LT.MAXX.AND.IY.LT.MAXX.AND.IZ.LT.MAXZ
     1             .AND.IX.GT.-MAXX.AND.IY.GT.-MAXX.AND.IZ.GT.-MAXZ)
     2             den=rho(ix,iy,iz)
          akg = 0.1727
          bkg = 0.333
         rnsg = den
         ecor = - akg*rnsg + (bkg*den)**2
         etotal = sqrt(etotal**2 + ecor)
          endif
          R(1,I) = R(1,I) + DT*P(1,I)/ETOTAL
          R(2,I) = R(2,I) + DT*P(2,I)/ETOTAL
          R(3,I) = R(3,I) + DT*P(3,I)/ETOTAL
            if ( cycbox.ne.0 ) then
              if ( r(1,i).gt. cycbox/2 ) r(1,i)=r(1,i)-cycbox
              if ( r(1,i).le.-cycbox/2 ) r(1,i)=r(1,i)+cycbox
              if ( r(2,i).gt. cycbox/2 ) r(2,i)=r(2,i)-cycbox
              if ( r(2,i).le.-cycbox/2 ) r(2,i)=r(2,i)+cycbox
              if ( r(3,i).gt. cycbox/2 ) r(3,i)=r(3,i)-cycbox
              if ( r(3,i).le.-cycbox/2 ) r(3,i)=r(3,i)+cycbox
            end if
          LB1=LB(I)
        IF(LB1.EQ.9)LD1=LD1+1
        IF(LB1.EQ.8)LD2=LD2+1
        IF(LB1.EQ.7)LD3=LD3+1
        IF(LB1.EQ.6)LD4=LD4+1
        IF(LB1.EQ.11)LN1=LN1+1
        IF(LB1.EQ.10)LN2=LN2+1
       IF((LB1.EQ.13).OR.(LB1.EQ.12))LN5=LN5+1
       IF(LB1.EQ.0)LE=LE+1
       IF(LB1.EQ.23)LKAON=LKAON+1
       IF(LB1.EQ.30)LKAONS=LKAONS+1
        IF(LB1.EQ.5)LP1=LP1+1
        IF(LB1.EQ.4)LP2=LP2+1
        IF(LB1.EQ.3)LP3=LP3+1
201     CONTINUE
        LP=LP1+LP2+LP3
        LD=LD1+LD2+LD3+LD4
        LN=LN1+LN2
        ALP=FLOAT(LP)/FLOAT(NUM)
        ALD=FLOAT(LD)/FLOAT(NUM)
        ALN=FLOAT(LN)/FLOAT(NUM)
       ALN5=FLOAT(LN5)/FLOAT(NUM)
        ATOTAL=ALP+ALD+ALN+0.5*ALN5
       ALE=FLOAT(LE)/FLOAT(NUM)
       ALKAON=FLOAT(LKAON)/FLOAT(NUM)
        if (icou .eq. 1) then
          iso=0
          do 1026 irun = 1,num
            iso=iso+massr(irun-1)
            do 1021 il = 1,massr(irun)
               temp(1,il) = 0.
               temp(2,il) = 0.
               temp(3,il) = 0.
 1021       continue
            do 1023 il = 1, massr(irun)
              i=iso+il
              if (zet(lb(i)).ne.0) then
                do 1022 jl = 1,il-1
                  j=iso+jl
                  if (zet(lb(j)).ne.0) then
                    ddx=r(1,i)-r(1,j)
                    ddy=r(2,i)-r(2,j)
                    ddz=r(3,i)-r(3,j)
                    rdiff = sqrt(ddx**2+ddy**2+ddz**2)
                    if (rdiff .le. 1.) rdiff = 1.
                    grp=zet(lb(i))*zet(lb(j))/rdiff**3
                    ddx=ddx*grp
                    ddy=ddy*grp
                    ddz=ddz*grp
                    temp(1,il)=temp(1,il)+ddx
                    temp(2,il)=temp(2,il)+ddy
                    temp(3,il)=temp(3,il)+ddz
                    temp(1,jl)=temp(1,jl)-ddx
                    temp(2,jl)=temp(2,jl)-ddy
                    temp(3,jl)=temp(3,jl)-ddz
                  end if
 1022          continue
              end if
 1023      continue
            do 1025 il = 1,massr(irun)
              i= iso+il
              if (zet(lb(i)).ne.0) then
                do 1024 idir = 1,3
                  p(idir,i) = p(idir,i) + temp(idir,il)
     &                                    * dt * 0.00144
 1024          continue
              end if
 1025      continue
 1026   continue
        end if
       spt=0
       spz=0
       ncen=0
       ekin=0
          NLOST = 0
          MEAN=0
         nquark=0
         nbaryn=0
           rads = 2.
           zras = 0.1
           denst = 0.
           edenst = 0.
          DO 6000 IRUN = 1,NUM
          MEAN=MEAN+MASSR(IRUN-1)
          DO 5800 J = 1,MASSR(irun)
          I=J+MEAN
           radut = sqrt(r(1,i)**2+r(2,i)**2)
       if( radut .le. rads )then
        if( abs(r(3,i)) .le. zras*nt*dt )then
         vols = 3.14159*rads**2*zras
         engs=sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+e(i)**2)
         gammas=1.
         if(e(i).ne.0.)gammas=engs/e(i)
         denst = denst + 1./gammas/vols
         edenst = edenst + engs/gammas/gammas/vols
        endif
       endif
         drr=sqrt(r(1,i)**2+r(2,i)**2+r(3,i)**2)
         if(drr.le.2.0)then
         spt=spt+p(1,i)**2+p(2,i)**2
         spz=spz+p(3,i)**2
         ncen=ncen+1
         ekin=ekin+sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+e(i)**2)-e(i)
         endif
              IX = NINT( R(1,I) )
              IY = NINT( R(2,I) )
              IZ = NINT( R(3,I) )
              IF(IX.LT.MAXX.AND.IY.LT.MAXX.AND.IZ.LT.MAXZ
     1          .AND.IX.GT.-MAXX.AND.IY.GT.-MAXX.AND.IZ.GT.-MAXZ) THEN
       if(rho(ix,iy,iz)/0.168.gt.dencut)go to 5800
       if((rho(ix,iy,iz)/0.168.gt.5.).and.(e(i).gt.0.9))
     &  nbaryn=nbaryn+1
       if(pel(ix,iy,iz).gt.2.0)nquark=nquark+1
       endif
        if(kpoten.ne.0.and.lb(i).eq.23)then
        den=0.
        IF(IX.LT.MAXX.AND.IY.LT.MAXX.AND.IZ.LT.MAXZ
     1       .AND.IX.GT.-MAXX.AND.IY.GT.-MAXX.AND.IZ.GT.-MAXZ) THEN
           den=rho(ix,iy,iz)
            akg = 0.1727
            bkg = 0.333
          rnsg = den
          ecor = - akg*rnsg + (bkg*den)**2
          etotal = sqrt(P(1,i)**2+p(2,I)**2+p(3,i)**2+e(i)**2 + ecor)
          ecor = - akg + 2.*bkg**2*den + 2.*bkg*etotal
        CALL GRADUK(IX,IY,IZ,GRADXk,GRADYk,GRADZk)
        P(1,I) = P(1,I) - DT * GRADXk*ecor/(2.*etotal)
        P(2,I) = P(2,I) - DT * GRADYk*ecor/(2.*etotal)
        P(3,I) = P(3,I) - DT * GRADZk*ecor/(2.*etotal)
        endif
         endif
        if(kpoten.ne.0.and.lb(i).eq.21)then
         den=0.
         IF(IX.LT.MAXX.AND.IY.LT.MAXX.AND.IZ.LT.MAXZ
     1        .AND.IX.GT.-MAXX.AND.IY.GT.-MAXX.AND.IZ.GT.-MAXZ) THEN
               den=rho(ix,iy,iz)
        CALL GRADUK(IX,IY,IZ,GRADXk,GRADYk,GRADZk)
            akg = 0.1727
            bkg = 0.333
          rnsg = den
          ecor = - akg*rnsg + (bkg*den)**2
          etotal = sqrt(P(1,i)**2+p(2,I)**2+p(3,i)**2+e(i)**2 + ecor)
          ecor = - akg + 2.*bkg**2*den - 2.*bkg*etotal
        P(1,I) = P(1,I) - DT * GRADXk*ecor/(2.*etotal)
        P(2,I) = P(2,I) - DT * GRADYk*ecor/(2.*etotal)
        P(3,I) = P(3,I) - DT * GRADZk*ecor/(2.*etotal)
        endif
         endif
       if(j.gt.mass)go to 5800         
        IF (ICOLL .NE. -1) THEN
           IF(IX.LT.MAXX.AND.IY.LT.MAXX.AND.IZ.LT.MAXZ
     1          .AND.IX.GT.-MAXX.AND.IY.GT.-MAXX.AND.IZ.GT.-MAXZ) THEN
                CALL GRADU(IPOT,IX,IY,IZ,GRADX,GRADY,GRADZ)
              TZ=0.
              GRADXN=0
              GRADYN=0
              GRADZN=0
              GRADXP=0
              GRADYP=0
              GRADZP=0
             IF(ICOU.EQ.1)THEN
                CALL GRADUP(IX,IY,IZ,GRADXP,GRADYP,GRADZP)
                CALL GRADUN(IX,IY,IZ,GRADXN,GRADYN,GRADZN)
               IF(ZET(LB(I)).NE.0)TZ=-1
               IF(ZET(LB(I)).EQ.0)TZ= 1
             END IF
           if(iabs(lb(i)).ge.14.and.iabs(lb(i)).le.17)then
              facl = 2./3.
            elseif(iabs(lb(i)).eq.40.or.iabs(lb(i)).eq.41)then
              facl = 1./3.
            else
              facl = 1.
            endif
        P(1,I) = P(1,I) - facl*DT * (GRADX+asy*(GRADXN-GRADXP)*TZ)
        P(2,I) = P(2,I) - facl*DT * (GRADY+asy*(GRADYN-GRADYP)*TZ)
        P(3,I) = P(3,I) - facl*DT * (GRADZ+asy*(GRADZN-GRADZP)*TZ)
                end if                                                       
              ENDIF
 5800       CONTINUE
 6000       CONTINUE
       CDEN=RHO(0,0,0)/0.168
        IF ((NT/NFREQ)*NFREQ .EQ. NT ) THEN
       if(icflow.eq.1)call flow(nt)
       endif
      IF (IAPAR2(1) .NE. 1) THEN
         CT = NT * DT
         IA = 0
         DO 1028 IRUN = 1, NUM
            DO 1027 IC = 1, MASSR(IRUN)
               IE = IA + IC
               RT(1, IC, IRUN) = R(1, IE)
               RT(2, IC, IRUN) = R(2, IE)
               RT(3, IC, IRUN) = R(3, IE)
               PT(1, IC, IRUN) = P(1, IE)
               PT(2, IC, IRUN) = P(2, IE)
               PT(3, IC, IRUN) = P(3, IE)
               ET(IC, IRUN) = E(IE)
               LT(IC, IRUN) = LB(IE)
               PROT(IC, IRUN) = PROPER(IE)
               dpertt(IC, IRUN)=dpertp(IE)
 1027       CONTINUE
            NP = MASSR(IRUN)
            NP1 = NPI(IRUN)
               ctlong = ct
             if(nt .eq. (ntmax-1))then
               ctlong = 1.E30
             elseif(nt .eq. ntmax)then
               go to 1111
             endif
            DO WHILE (NP1.LE.MULTI1(IRUN).AND.
     &           FT1(NP1, IRUN) .GT. ((NT-1) * DT) .AND. 
     &           FT1(NP1, IRUN) .LE. ctlong)
               NP = NP + 1
               UDT = (CT - FT1(NP1, IRUN)) / EE1(NP1, IRUN)
               if(nt.eq.(ntmax-1)) then
                  ftsvt(NP,IRUN)=FT1(NP1, IRUN)
                  if(FT1(NP1, IRUN).gt.ct) UDT=0.
               endif
               RT(1, NP, IRUN) = GX1(NP1, IRUN) + 
     &              PX1(NP1, IRUN) * UDT
               RT(2, NP, IRUN) = GY1(NP1, IRUN) + 
     &              PY1(NP1, IRUN) * UDT
               RT(3, NP, IRUN) = GZ1(NP1, IRUN) + 
     &              PZ1(NP1, IRUN) * UDT
               PT(1, NP, IRUN) = PX1(NP1, IRUN)
               PT(2, NP, IRUN) = PY1(NP1, IRUN)
               PT(3, NP, IRUN) = PZ1(NP1, IRUN)
               ET(NP, IRUN) = XM1(NP1, IRUN)
               LT(NP, IRUN) = IARFLV(ITYP1(NP1, IRUN))
               dpertt(NP,IRUN)=dpp1(NP1,IRUN)
               NP1 = NP1 + 1
               PROT(NP, IRUN) = 1.
            END DO
 1111      continue
            NPI(IRUN) = NP1
            IA = IA + MASSR(IRUN)
            MASSR(IRUN) = NP
 1028    CONTINUE
         IA = 0
         DO 1030 IRUN = 1, NUM
            IA = IA + MASSR(IRUN - 1)
            DO 1029 IC = 1, MASSR(IRUN)
               IE = IA + IC
               R(1, IE) = RT(1, IC, IRUN)
               R(2, IE) = RT(2, IC, IRUN)
               R(3, IE) = RT(3, IC, IRUN)
               P(1, IE) = PT(1, IC, IRUN)
               P(2, IE) = PT(2, IC, IRUN)
               P(3, IE) = PT(3, IC, IRUN)
               E(IE) = ET(IC, IRUN)
               LB(IE) = LT(IC, IRUN)
               PROPER(IE) = PROT(IC, IRUN)
               if(nt.eq.(ntmax-1)) ftsv(IE)=ftsvt(IC,IRUN)
               dpertp(IE)=dpertt(IC, IRUN)
 1029       CONTINUE
            call hbtout(MASSR(IRUN),nt,ntmax)
 1030    CONTINUE
      END IF
10000       continue
        iss=0
        do 1032 lrun=1,num
           iss=iss+massr(lrun-1)
           do 1031 l0=1,massr(lrun)
              ipart=iss+l0
 1031      continue
 1032   continue
      IF (IAPAR2(1) .NE. 1) THEN
        IA = 0
        DO 1035 IRUN = 1, NUM
           IA = IA + MASSR(IRUN - 1)
           NP1 = NPI(IRUN)
           NSH = MASSR(IRUN) - NP1 + 1
           MULTI1(IRUN) = MULTI1(IRUN) + NSH
           IF (NSH .GT. 0) THEN
              IB = MULTI1(IRUN)
              IE = MASSR(IRUN) + 1
              II = -1
           ELSE IF (NSH .LT. 0) THEN
              IB = MASSR(IRUN) + 1
              IE = MULTI1(IRUN)
              II = 1
           END IF
           IF (NSH .NE. 0) THEN
              DO 1033 I = IB, IE, II
                 J = I - NSH
                 ITYP1(I, IRUN) = ITYP1(J, IRUN)
                 GX1(I, IRUN) = GX1(J, IRUN)
                 GY1(I, IRUN) = GY1(J, IRUN)
                 GZ1(I, IRUN) = GZ1(J, IRUN)
                 FT1(I, IRUN) = FT1(J, IRUN)
                 PX1(I, IRUN) = PX1(J, IRUN)
                 PY1(I, IRUN) = PY1(J, IRUN)
                 PZ1(I, IRUN) = PZ1(J, IRUN)
                 EE1(I, IRUN) = EE1(J, IRUN)
                 XM1(I, IRUN) = XM1(J, IRUN)
                 PRO1(I, IRUN) = PRO1(J, IRUN)
                 dpp1(I,IRUN)=dpp1(J,IRUN)
 1033         CONTINUE
           END IF
           DO 1034 I = 1, MASSR(IRUN)
              IB = IA + I
              ITYP1(I, IRUN) = INVFLV(LB(IB))
              GX1(I, IRUN) = R(1, IB)
              GY1(I, IRUN) = R(2, IB)
              GZ1(I, IRUN) = R(3, IB)
              if(FT1(I, IRUN).lt.CT) FT1(I, IRUN) = CT
              PX1(I, IRUN) = P(1, IB)
              PY1(I, IRUN) = P(2, IB)
              PZ1(I, IRUN) = P(3, IB)
              XM1(I, IRUN) = E(IB)
              EE1(I, IRUN) = SQRT(PX1(I, IRUN) ** 2 + 
     &             PY1(I, IRUN) ** 2 +
     &             PZ1(I, IRUN) ** 2 + 
     &             XM1(I, IRUN) ** 2)
              PRO1(I, IRUN) = PROPER(IB)
 1034      CONTINUE
 1035   CONTINUE
      END IF
50000   CONTINUE
      RETURN
      END
