      SUBROUTINE RELCOL(LCOLL,LBLOC,LCNNE,LDD,LPP,lppk,
     &LPN,lpd,lrho,lomega,LKN,LNNK,LDDK,LNDK,LCNND,LCNDN,
     &LDIRT,LDECAY,LRES,LDOU,LDDRHO,LNNRHO,LNNOM,
     &NT,ntmax,sp,akaon,sk)
      PARAMETER      (MAXSTR=150001,MAXR=1,PI=3.1415926)
      parameter      (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
      PARAMETER      (AKA=0.498,ALA=1.1157,ASA=1.1974,aks=0.895)
      PARAMETER      (AA1=1.26,APHI=1.02,AP1=0.13496)
      parameter            (maxx=20,maxz=24)
      parameter            (rrkk=0.6,prkk=0.3,srhoks=5.,ESBIN=0.04)
      DIMENSION MASSRN(0:MAXR),RT(3,MAXSTR),PT(3,MAXSTR),ET(MAXSTR)
      DIMENSION LT(MAXSTR), PROT(MAXSTR)
      COMMON   /AA/  R(3,MAXSTR)
      COMMON   /BB/  P(3,MAXSTR)
      COMMON   /CC/  E(MAXSTR)
      COMMON  /DD/      RHO(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHOP(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHON(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ)
      COMMON   /EE/  ID(MAXSTR),LB(MAXSTR)
      COMMON   /HH/  PROPER(MAXSTR)
      common /ff/f(-mx:mx,-my:my,-mz:mz,-mpx:mpx,-mpy:mpy,-mpz:mpzp)
      common   /gg/  dx,dy,dz,dpx,dpy,dpz
      COMMON   /INPUT/ NSTAR,NDIRCT,DIR
      COMMON   /NN/NNN
      COMMON   /RR/  MASSR(0:MAXR)
      common   /ss/  inout(20)
      COMMON   /BG/BETAX,BETAY,BETAZ,GAMMA
      COMMON   /RUN/NUM
      COMMON   /PA/RPION(3,MAXSTR,MAXR)
      COMMON   /PB/PPION(3,MAXSTR,MAXR)
      COMMON   /PC/EPION(MAXSTR,MAXR)
      COMMON   /PD/LPION(MAXSTR,MAXR)
      COMMON   /PE/PROPI(MAXSTR,MAXR)
      COMMON   /KKK/TKAON(7),EKAON(7,0:2000)
      COMMON  /KAON/    AK(3,50,36),SPECK(50,36,7),MF
      COMMON/TABLE/ xarray(0:1000),earray(0:1000)
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      common/leadng/lb1,px1,py1,pz1,em1,e1,xfnl,yfnl,zfnl,tfnl,
     1 px1n,py1n,pz1n,dp1n
      COMMON/tdecay/tfdcy(MAXSTR),tfdpi(MAXSTR,MAXR),tft(MAXSTR)
      common /lastt/itimeh,bimp 
      COMMON/ppbmas/niso(15),nstate,ppbm(15,2),thresh(15),weight(15)
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
      COMMON/hbt/lblast(MAXSTR),xlast(4,MAXSTR),plast(4,MAXSTR),nlast
      common/resdcy/NSAV,iksdcy
      COMMON/RNDF77/NSEED
      COMMON/FTMAX/ftsv(MAXSTR),ftsvt(MAXSTR, MAXR)
      dimension ftpisv(MAXSTR,MAXR),fttemp(MAXSTR)
      common /dpi/em2,lb2
      common/phidcy/iphidcy,pttrig,ntrig,maxmiss,ipi0dcy
      DIMENSION dptemp(MAXSTR)
      common /para8/ idpert,npertd,idxsec
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
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
      call inidcy
      RESONA=5.
      NODELT=0
      SUMSRT =0.
      LCOLL  = 0
      LBLOC  = 0
      LCNNE  = 0
      LDD  = 0
      LPP  = 0
      lpd  = 0
      lpdr=0
      lrho = 0
      lrhor=0
      lomega=0
      lomgar=0
      LPN  = 0
      LKN  = 0
      LNNK = 0
      LDDK = 0
      LNDK = 0
      lppk =0
      LCNND  = 0
      LCNDN  = 0
      LDIRT  = 0
      LDECAY = 0
      LRES   = 0
      Ldou   = 0
      LDDRHO = 0
      LNNRHO = 0
      LNNOM  = 0
      MSUM   = 0
      MASSRN(0)=0
      DO 1002 IL=1,5
         TKAON(IL)=0
         DO 1001 IS=1,2000
            EKAON(IL,IS)=0
 1001    CONTINUE
 1002 CONTINUE
      DO 1004 i =1,NUM
         DO 1003 j =1,MAXSTR
            PROPI(j,i) = 1.
 1003    CONTINUE
 1004 CONTINUE
      do 1102 i=1,maxstr
         fttemp(i)=0.
         do 1101 irun=1,maxr
            ftpisv(i,irun)=0.
 1101    continue
 1102 continue
      sp=0
      akaon=0
      sk=0
      MASS = 0
      DO 1000 IRUN = 1,NUM
         NNN=0
         MSUM=MSUM+MASSR(IRUN-1)
         J10=2
         IF(NT.EQ.NTMAX)J10=1
         DO 800 J1 = J10,MASSR(IRUN)
            I1  = J1 + MSUM
            IF(E(I1).EQ.0.)GO TO 798
            IF (LB(I1) .LT. -45 .OR. LB(I1) .GT. 45) GOTO 798
            X1  = R(1,I1)
            Y1  = R(2,I1)
            Z1  = R(3,I1)
            PX1 = P(1,I1)
            PY1 = P(2,I1)
            PZ1 = P(3,I1)
            EM1 = E(I1)
            am1= em1
            E1  = SQRT( EM1**2 + PX1**2 + PY1**2 + PZ1**2 )
            ID1 = ID(I1)
            LB1 = LB(I1)
            if(nt.eq.ntmax.and.(lb1.eq.21.or.lb1.eq.23)) then
               pk0=RANART(NSEED)
               if(pk0.lt.0.25) then
                  LB(I1)=22
               elseif(pk0.lt.0.50) then
                  LB(I1)=24
               endif
               LB1=LB(I1)
            endif
            if(lb1.eq.0.or.lb1.eq.25.or.lb1.eq.26.or.lb1.eq.27
     &           .or.lb1.eq.28.or.lb1.eq.29.or.iabs(lb1).eq.30
     &           .or.(iabs(lb1).ge.6.and.iabs(lb1).le.13)
     &           .or.(iksdcy.eq.1.and.lb1.eq.24)
     &           .or.iabs(lb1).eq.16
     &           .or.(ipi0dcy.eq.1.and.nt.eq.ntmax.and.lb1.eq.4)) then
               continue
            else
               goto 1
            endif
         IF(lb1.ge.25.and.lb1.le.27) then
             wid=0.151
         ELSEIF(lb1.eq.28) then
             wid=0.00841
         ELSEIF(lb1.eq.29) then
             wid=0.00443
          ELSEIF(iabs(LB1).eq.30) then
             WID=0.051
         ELSEIF(lb1.eq.0) then
             wid=1.18e-6
         ELSEIF(iksdcy.eq.1.and.lb1.eq.24) then
             wid=7.36e-15
         ELSEIF(iabs(lb1).eq.16) then
             wid=8.87e-6
          ELSEIF(LB1.EQ.32) then
             call WIDA1(EM1,rhomp,WID,iseed)
          ELSEIF(iabs(LB1).ge.6.and.iabs(LB1).le.9) then
             WID=WIDTH(EM1)
          ELSEIF((iabs(LB1).EQ.10).OR.(iabs(LB1).EQ.11)) then
             WID=W1440(EM1)
          ELSEIF((iabs(LB1).EQ.12).OR.(iabs(LB1).EQ.13)) then
             WID=W1535(EM1)
          ELSEIF(ipi0dcy.eq.1.and.nt.eq.ntmax.and.lb1.eq.4) then
             wid=7.85e-9
          ENDIF
          if(nt.eq.ntmax)then
             pdecay=1.1
             if(iphidcy.eq.0.and.iabs(LB1).eq.29) pdecay=0.
          else
             T0=0.19733/WID
             GFACTR=E1/EM1
             T0=T0*GFACTR
             IF(T0.GT.0.)THEN
                PDECAY=1.-EXP(-DT/T0)
             ELSE
                PDECAY=0.
             ENDIF
          endif
          XDECAY=RANART(NSEED)
          IF(XDECAY.LT.PDECAY) THEN
             idecay=irun
             tfnl=nt*dt
             if(nt.eq.ntmax.and.ftsv(i1).gt.((ntmax-1)*dt)) 
     1            tfnl=ftsv(i1)
             xfnl=x1
             yfnl=y1
             zfnl=z1
             if(lb1.eq.0.or.lb1.eq.25.or.lb1.eq.26.or.lb1.eq.27
     &           .or.lb1.eq.28.or.lb1.eq.29.or.iabs(lb1).eq.30
     &           .or.(iabs(lb1).ge.6.and.iabs(lb1).le.9)
     &           .or.(iksdcy.eq.1.and.lb1.eq.24)
     &           .or.iabs(lb1).eq.16
     &           .or.(ipi0dcy.eq.1.and.nt.eq.ntmax.and.lb1.eq.4)) then
                call resdec(i1,nt,nnn,wid,idecay,0)
                p(1,i1)=px1n
                p(2,i1)=py1n
                p(3,i1)=pz1n
                dpertp(i1)=dp1n
                if(nt.eq.ntmax) then
                   R(1,i1)=xfnl
                   R(2,i1)=yfnl
                   R(3,i1)=zfnl
                   tfdcy(i1)=tfnl
                endif
                if(iabs(lb1).ge.6.and.iabs(lb1).le.9) then
                   LDECAY=LDECAY+1
                endif
             elseif(iabs(LB1).EQ.10.OR.iabs(LB1).EQ.11) THEN
                NNN=NNN+1
                LDECAY=LDECAY+1
                PNSTAR=1.
                IF(E(I1).GT.1.22)PNSTAR=0.6
                IF(RANART(NSEED).LE.PNSTAR)THEN
                   CALL DECAY(idecay,I1,NNN,ISEED,wid,nt)
                ELSE
                   CALL DECAY2(idecay,I1,NNN,ISEED,wid,nt)
                   NNN=NNN+1
                ENDIF
             elseif(iabs(LB1).eq.12.or.iabs(LB1).eq.13) then
                NNN=NNN+1
                CALL DECAY(idecay,I1,NNN,ISEED,wid,nt)
                LDECAY=LDECAY+1
             endif
             if(nt.eq.ntmax) then
                if(lb(i1).eq.25.or.lb(i1).eq.26.or.lb(i1).eq.27) then
                   wid=0.151
                elseif(lb(i1).eq.0) then
                   wid=1.18e-6
                elseif(lb(i1).eq.24.and.iksdcy.eq.1) then
                   wid=7.36e-15
                elseif(ipi0dcy.eq.1.and.lb(i1).eq.4) then
                   wid=7.85e-9
                else
                   goto 9000
                endif
                LB1=LB(I1)
                PX1=P(1,I1)
                PY1=P(2,I1)
                PZ1=P(3,I1)
                EM1=E(I1)
                E1=SQRT(EM1**2+PX1**2+PY1**2+PZ1**2)
                call resdec(i1,nt,nnn,wid,idecay,0)
                p(1,i1)=px1n
                p(2,i1)=py1n
                p(3,i1)=pz1n
                R(1,i1)=xfnl
                R(2,i1)=yfnl
                R(3,i1)=zfnl
                tfdcy(i1)=tfnl
                dpertp(i1)=dp1n
             endif
             if(nt.eq.ntmax.and.ipi0dcy.eq.1.and.lb(i1).eq.4) then
                wid=7.85e-9
                LB1=LB(I1)
                PX1=P(1,I1)
                PY1=P(2,I1)
                PZ1=P(3,I1)
                EM1=E(I1)
                E1=SQRT(EM1**2+PX1**2+PY1**2+PZ1**2)
                call resdec(i1,nt,nnn,wid,idecay,0)
                p(1,i1)=px1n
                p(2,i1)=py1n
                p(3,i1)=pz1n
                R(1,i1)=xfnl
                R(2,i1)=yfnl
                R(3,i1)=zfnl
                tfdcy(i1)=tfnl
                dpertp(i1)=dp1n
             endif
 9000        go to 798
          ENDIF
 1        if(nt.eq.ntmax)go to 798
          X1 = R(1,I1)
          Y1 = R(2,I1)
          Z1 = R(3,I1)
           DO 600 J2 = 1,J1-1
            I2  = J2 + MSUM
            IF(E(I2).EQ.0.) GO TO 600
            IF(E(I1).EQ.0.) GO TO 800
            IF (LB(I2) .LT. -45 .OR. LB(I2) .GT. 45) GOTO 600
            X2=R(1,I2)
            Y2=R(2,I2)
            Z2=R(3,I2)
            dr0max=5.
            ilb1=iabs(LB(I1))
            ilb2=iabs(LB(I2))
            IF(ilb1.EQ.42.or.ilb2.EQ.42) THEN
               if((ILB1.GE.1.AND.ILB1.LE.2)
     1              .or.(ILB1.GE.6.AND.ILB1.LE.13)
     2              .or.(ILB2.GE.1.AND.ILB2.LE.2)
     3              .or.(ILB2.GE.6.AND.ILB2.LE.13)) then
                  if((lb(i1)*lb(i2)).gt.0) dr0max=10.
               endif
            ENDIF
            if(((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2).GT.dr0max**2)
     1           GO TO 600
            IF (ID(I1)*ID(I2).EQ.IAVOID) GOTO 400
            ID1=ID(I1)
            ID2 = ID(I2)
            ix1= nint(x1/dx)
            iy1= nint(y1/dy)
            iz1= nint(z1/dz)
            PX1=P(1,I1)
            PY1=P(2,I1)
            PZ1=P(3,I1)
            EM1=E(I1)
            AM1=EM1
            LB1=LB(I1)
            E1=SQRT(EM1**2+PX1**2+PY1**2+PZ1**2)
            IPX1=NINT(PX1/DPX)
            IPY1=NINT(PY1/DPY)
            IPZ1=NINT(PZ1/DPZ)         
            LB2 = LB(I2)
            PX2 = P(1,I2)
            PY2 = P(2,I2)
            PZ2 = P(3,I2)
            EM2=E(I2)
            AM2=EM2
            lb1i=lb(i1)
            lb2i=lb(i2)
            px1i=P(1,I1)
            py1i=P(2,I1)
            pz1i=P(3,I1)
            em1i=E(I1)
            px2i=P(1,I2)
            py2i=P(2,I2)
            pz2i=P(3,I2)
            em2i=E(I2)
            eini=SQRT(E(I1)**2+P(1,I1)**2+P(2,I1)**2+P(3,I1)**2)
     1           +SQRT(E(I2)**2+P(1,I2)**2+P(2,I2)**2+P(3,I2)**2)
            pxini=P(1,I1)+P(1,I2)
            pyini=P(2,I1)+P(2,I2)
            pzini=P(3,I1)+P(3,I2)
            nnnini=nnn
            iblock=0
            DELTR0=3.
        if( (iabs(lb1).ge.14.and.iabs(lb1).le.17) .or.
     &      (iabs(lb1).ge.30.and.iabs(lb1).le.45) ) DELTR0=5.0
        if( (iabs(lb2).ge.14.and.iabs(lb2).le.17) .or.
     &      (iabs(lb2).ge.30.and.iabs(lb2).le.45) ) DELTR0=5.0
            if(lb1.eq.28.and.lb2.eq.28) DELTR0=4.84
            if((lb1.ge.3.and.lb1.le.5).and.(lb2.ge.3.and.lb2.le.5)) then
               E2=SQRT(EM2**2+PX2**2+PY2**2+PZ2**2)
         spipi=(e1+e2)**2-(px1+px2)**2-(py1+py2)**2-(pz1+pz2)**2
               if(spipi.ge.(4*0.77**2)) DELTR0=3.5
            endif
        IF (LB1.EQ.23 .AND. (LB2.GE.14.AND.LB2.LE.17)) GOTO 3699
        IF (LB2.EQ.23 .AND. (LB1.GE.14.AND.LB1.LE.17)) GOTO 3699
       if(lb1.eq.21.and.lb2.eq.23)go to 3699
       if(lb2.eq.21.and.lb1.eq.23)go to 3699
       if(lb1.eq.30.and.lb2.eq.21)go to 3699
       if(lb2.eq.30.and.lb1.eq.21)go to 3699
       if(lb1.eq.-30.and.lb2.eq.23)go to 3699
       if(lb2.eq.-30.and.lb1.eq.23)go to 3699
       if(lb1.eq.-30.and.lb2.eq.30)go to 3699
       if(lb2.eq.-30.and.lb1.eq.30)go to 3699
      if(lb1.eq.21.or.lb1.eq.23) then
         if(lb2.eq.0.or.(lb2.ge.25.and.lb2.le.28)) then
            go to 3699
         endif
      elseif(lb2.eq.21.or.lb2.eq.23) then
         if(lb1.eq.0.or.(lb1.ge.25.and.lb1.le.28)) then
            goto 3699
         endif
      endif
      if(iabs(lb1).eq.30 .and.
     1     (lb2.eq.0.or.(lb2.ge.25.and.lb2.le.28)
     2     .or.(lb2.ge.3.and.lb2.le.5))) then
         go to 3699
      elseif(iabs(lb2).eq.30 .and.
     1        (lb1.eq.0.or.(lb1.ge.25.and.lb1.le.28)
     2        .or.(lb1.ge.3.and.lb1.le.5))) then
         goto 3699
        elseif( iabs(lb1).eq.30 .and.
     1     (iabs(lb2).eq.1.or.iabs(lb2).eq.2.or.
     2     (iabs(lb2).ge.6.and.iabs(lb2).le.13)) )then
              go to 3699
           endif
         if( iabs(lb2).eq.30 .and.
     1         (iabs(lb1).eq.1.or.iabs(lb1).eq.2.or.
     2         (iabs(lb1).ge.6.and.iabs(lb1).le.13)) )then
                go to 3699
        endif                                                              
        if((lb1.eq.23.or.lb1.eq.21).and.
     1       (iabs(lb2).eq.1.or.iabs(lb2).eq.2.or.
     2       (iabs(lb2).ge.6.and.iabs(lb2).le.13))) then
           go to 3699
        elseif((lb2.eq.23.or.lb2.eq.21).and.
     1       (iabs(lb1).eq.1.or.iabs(lb1).eq.2.or.
     2       (iabs(lb1).ge.6.and.iabs(lb1).le.13))) then
           go to 3699
        endif
        rppmax=3.57   
        if((lb1.eq.-1.or.lb1.eq.-2.or.(lb1.ge.-13.and.lb1.le.-6))
     1 .and.(lb2.eq.1.or.lb2.eq.2.or.(lb2.ge.6.and.lb2.le.13))) then
            DELTR0 = RPPMAX
            GOTO 2699
       else if((lb2.eq.-1.or.lb2.eq.-2.or.(lb2.ge.-13.and.lb2.le.-6))
     1 .and.(lb1.eq.1.or.lb1.eq.2.or.(lb1.ge.6.and.lb1.le.13))) then
            DELTR0 = RPPMAX
            GOTO 2699
         END IF
        if( (iabs(lb1).ge.14.and.iabs(lb1).le.17) .or.
     &      (iabs(lb2).ge.14.and.iabs(lb2).le.17) )go to 3699
         IF (iabs(LB1).EQ.42.or.iabs(LB2).EQ.42) THEN
            ilb1=iabs(LB1)
            ilb2=iabs(LB2)
            if((ILB1.GE.1.AND.ILB1.LE.2)
     1           .or.(ILB1.GE.6.AND.ILB1.LE.13)
     2           .or.(ILB2.GE.1.AND.ILB2.LE.2)
     3           .or.(ILB2.GE.6.AND.ILB2.LE.13)) then
               if((lb1*lb2).gt.0) deltr0=9.5
            endif
         ENDIF
        if( (iabs(lb1).ge.40.and.iabs(lb1).le.45) .or. 
     &      (iabs(lb2).ge.40.and.iabs(lb2).le.45) )go to 3699
         IF( (lb1.eq.29 .and.((lb2.ge.1.and.lb2.le.13).or.  
     &       (lb2.ge.21.and.lb2.le.28).or.iabs(lb2).eq.30)) .OR.
     &     (lb2.eq.29 .and.((lb1.ge.1.and.lb1.le.13).or.
     &       (lb1.ge.21.and.lb1.le.28).or.iabs(lb1).eq.30)) )THEN
             DELTR0=3.0
             go to 3699
        endif
        If(iabs(lb1).eq.30.or.iabs(lb2).eq.30) go to 400
         If(lb1.eq.23.and.(lb2.lt.1.or.lb2.gt.17))go to 400
         If(lb2.eq.23.and.(lb1.lt.1.or.lb1.gt.17))go to 400
         if( ((lb1.le.-1.and.lb1.ge.-13)
     &        .and.(lb2.eq.0.or.(lb2.ge.3.and.lb2.le.5)
     &            .or.(lb2.ge.25.and.lb2.le.28))) 
     &      .OR.((lb2.le.-1.and.lb2.ge.-13)
     &         .and.(lb1.eq.0.or.(lb1.ge.3.and.lb1.le.5)
     &              .or.(lb1.ge.25.and.lb1.le.28))) ) then
         elseIF( ((LB1.eq.-1.or.lb1.eq.-2).
     &             and.(LB2.LT.-5.and.lb2.ge.-13))
     &      .OR. ((LB2.eq.-1.or.lb2.eq.-2).
     &             and.(LB1.LT.-5.and.lb1.ge.-13)) )then
         elseIF((LB1.eq.-1.or.lb1.eq.-2)
     &     .AND.(LB2.eq.-1.or.lb2.eq.-2))then
         elseIF((LB1.LT.-5.and.lb1.ge.-13).AND.
     &          (LB2.LT.-5.and.lb2.ge.-13)) then
       endif               
 2699    CONTINUE
         IF (LB1 .EQ. 1 .OR. LB1 .EQ. 2 .OR. (LB1 .GE. 6 .AND.
     &        LB1 .LE. 17)) THEN
            IF (LB2 .EQ. 1 .OR. LB2 .EQ. 2 .OR. (LB2 .GE. 6 .AND.
     &           LB2 .LE. 17)) THEN
               DELTR0 = 2.
            END IF
         END IF
 3699   RSQARE = (X1-X2)**2 + (Y1-Y2)**2 + (Z1-Z2)**2
        IF (RSQARE .GT. DELTR0**2) GO TO 400
            ix2 = nint(x2/dx)
            iy2 = nint(y2/dy)
            iz2 = nint(z2/dz)
            ipx2 = nint(px2/dpx)
            ipy2 = nint(py2/dpy)
            ipz2 = nint(pz2/dpz)
            CALL CMS(I1,I2,PCX,PCY,PCZ,SRT)
          drmax=dr0max
          call distc0(drmax,deltr0,DT,
     1         Ifirst,PCX,PCY,PCZ,
     2         x1,y1,z1,px1,py1,pz1,em1,x2,y2,z2,px2,py2,pz2,em2)
          if(Ifirst.eq.-1) goto 400
         ISS=NINT(SRT/ESBIN)
         if(ISS.gt.2000) ISS=2000
         IF (iabs(LB1).EQ.42.or.iabs(LB2).EQ.42) THEN
            ilb1=iabs(LB1)
            ilb2=iabs(LB2)
            if(LB1.eq.0.or.(LB1.GE.3.AND.LB1.LE.5)
     1           .or.(LB1.GE.25.AND.LB1.LE.28)
     2           .or.
     3           LB2.eq.0.or.(LB2.GE.3.AND.LB2.LE.5)
     4           .or.(LB2.GE.25.AND.LB2.LE.28)) then
               GOTO 505
            elseif(((ILB1.GE.1.AND.ILB1.LE.2)
     1              .or.(ILB1.GE.6.AND.ILB1.LE.13)
     2              .or.(ILB2.GE.1.AND.ILB2.LE.2)
     3              .or.(ILB2.GE.6.AND.ILB2.LE.13))
     4              .and.(lb1*lb2).gt.0) then
               GOTO 506
            else
               GOTO 400
            endif
         ENDIF
          if( ((lb1.eq.23.or.lb1.eq.30).and.
     &         (lb2.eq.-1.or.lb2.eq.-2.or.(lb2.ge.-13.and.lb2.le.-6))) 
     &         .OR.((lb2.eq.23.or.lb2.eq.30).and.
     &         (lb1.eq.-1.or.lb1.eq.-2.or.(lb1.ge.-13.and.lb1.le.-6))) )
     &         then
             bmass=0.938
             if(srt.le.(bmass+aka)) then
                pkaon=0.
             else
                pkaon=sqrt(((srt**2-(aka**2+bmass**2))
     1               /2./bmass)**2-aka**2)
             endif
             sigela = 0.5 * (AKPEL(PKAON) + AKNEL(PKAON))
             SIGSGM = 1.5 * AKPSGM(PKAON) + AKNSGM(PKAON)
             SIG = sigela + SIGSGM + AKPLAM(PKAON)
             if(sig.gt.1.e-7) then
                icase=3
                brel=sigela/sig
                brsgm=sigsgm/sig
                brsig = sig
                nchrg = 1
                go to 3555
             endif
             go to 400
          endif
          if(((lb1.ge.-17.and.lb1.le.-14).and.(lb2.ge.3.and.lb2.le.5)) 
     &         .OR.((lb2.ge.-17.and.lb2.le.-14)
     &         .and.(lb1.ge.3.and.lb1.le.5)))then
             nchrg=-100
             if((lb1.eq.-15.and.(lb2.eq.5.or.lb2.eq.27)).OR.
     &            (lb2.eq.-15.and.(lb1.eq.5.or.lb1.eq.27))) then
                nchrg=-2
                bmass=1.232
                go to 110
             endif
             if( (lb1.eq.-15.and.(lb2.eq.0.or.lb2.eq.4.or.lb2.eq.26.or.
     &            lb2.eq.28)).OR.(lb2.eq.-15.and.(lb1.eq.0.or.
     &            lb1.eq.4.or.lb1.eq.26.or.lb1.eq.28)).OR.
     &   ((lb1.eq.-14.or.lb1.eq.-16).and.(lb2.eq.5.or.lb2.eq.27)).OR.
     &   ((lb2.eq.-14.or.lb2.eq.-16).and.(lb1.eq.5.or.lb1.eq.27)) )then
                nchrg=-1
                bmass=0.938
                go to 110
             endif
             if(  (lb1.eq.-15.and.(lb2.eq.3.or.lb2.eq.25)).OR.
     &            (lb2.eq.-15.and.(lb1.eq.3.or.lb1.eq.25)).OR.
     &            (lb1.eq.-17.and.(lb2.eq.5.or.lb2.eq.27)).OR.
     &            (lb2.eq.-17.and.(lb1.eq.5.or.lb1.eq.27)).OR.
     &            ((lb1.eq.-14.or.lb1.eq.-16).and.(lb2.eq.0.or.lb2.eq.4
     &            .or.lb2.eq.26.or.lb2.eq.28)).OR.
     &            ((lb2.eq.-14.or.lb2.eq.-16).and.(lb1.eq.0.or.lb1.eq.4
     &            .or.lb1.eq.26.or.lb1.eq.28)) )then
               nchrg=0
                bmass=0.938
                go to 110
             endif
             if( (lb1.eq.-17.and.(lb2.eq.0.or.lb2.eq.4.or.lb2.eq.26.or.
     &            lb2.eq.28)).OR.(lb2.eq.-17.and.(lb1.eq.0.or.
     &            lb1.eq.4.or.lb1.eq.26.or.lb1.eq.28)).OR.
     &  ((lb1.eq.-14.or.lb1.eq.-16).and.(lb2.eq.3.or.lb2.eq.25)).OR.
     &  ((lb2.eq.-14.or.lb2.eq.-16).and.(lb1.eq.3.or.lb1.eq.25)))then
               nchrg=1
                bmass=1.232
             endif
 110         sig = 0.
         if(nchrg.ne.-100.and.srt.ge.(aka+bmass))then
            icase=4
            pkaon=sqrt(((srt**2-(aka**2+0.938**2))/2./0.938)**2-aka**2)
            if(lb1.eq.-14.or.lb2.eq.-14) then
               if(nchrg.ge.0) sigma0=akPlam(pkaon)
               if(nchrg.lt.0) sigma0=akNlam(pkaon)
            else
               if(nchrg.ge.0) sigma0=akPsgm(pkaon)
               if(nchrg.lt.0) sigma0=akNsgm(pkaon)
               SIGMA0 = 1.5 * AKPSGM(PKAON) + AKNSGM(PKAON)
            endif
            sig=(srt**2-(aka+bmass)**2)*(srt**2-(aka-bmass)**2)/
     &           (srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)*sigma0
            if(nchrg.eq.-2.or.nchrg.eq.2) sig=2.*sig
            IF (LB1 .EQ. -14 .OR. LB2 .EQ. -14) THEN
               SIG = 4.0 / 3.0 * SIG
            ELSE IF (NCHRG .EQ. -2 .OR. NCHRG .EQ. 2) THEN
               SIG = 8.0 / 9.0 * SIG
            ELSE
               SIG = 4.0 / 9.0 * SIG
            END IF
         endif
         icase=4
         sigela = 10.
         sig = sig + sigela
         brel= sigela/sig
         brsgm=0.
         brsig = sig
         go to 3555
      endif
      if( ((lb1.eq.21.or.lb1.eq.-30).and.(lb2.ge.14.and.lb2.le.17)).OR.
     &  ((lb2.eq.21.or.lb2.eq.-30).and.(lb1.ge.14.and.lb1.le.17)) )then
          kp = 0
          go to 3455
        endif
      if( ((lb1.eq.23.or.lb1.eq.30).and.(lb2.le.-14.and.lb2.ge.-17)).OR.
     &  ((lb2.eq.23.or.lb2.eq.30).and.(lb1.le.-14.and.lb1.ge.-17)) )then
          kp = 1
          go to 3455
        endif
       if( ((lb1.eq.21.or.lb1.eq.-30).and.(lb2.eq.40.or.lb2.eq.41)).OR.
     & ((lb2.eq.21.or.lb2.eq.-30).and.(lb1.eq.40.or.lb1.eq.41)) )then
          kp = 0
          go to 3455
        endif
       if( ((lb1.eq.23.or.lb1.eq.30).and.(lb2.eq.-40.or.lb2.eq.-41)).OR.
     &  ((lb2.eq.23.or.lb2.eq.30).and.(lb1.eq.-40.or.lb1.eq.-41)) )then
          kp = 1
          go to 3455
        endif
         kp = 3
       if( (((lb1.ge.3.and.lb1.le.5).or.lb1.eq.0) 
     &       .and.(iabs(lb2).eq.40.or.iabs(lb2).eq.41))
     & .OR. (((lb2.ge.3.and.lb2.le.5).or.lb2.eq.0) 
     &       .and.(iabs(lb1).eq.40.or.iabs(lb1).eq.41)) )go to 3455
        if( ((lb1.ge.3.and.lb1.le.5).and.iabs(lb2).eq.45)
     &  .OR.((lb2.ge.3.and.lb2.le.5).and.iabs(lb1).eq.45) )go to 3455
        IF (LB1.EQ.23 .AND. (LB2.GE.14.AND.LB2.LE.17)) GOTO 5699
        IF (LB2.EQ.23 .AND. (LB1.GE.14.AND.LB1.LE.17)) GOTO 5699
       IF (LB1.EQ.21 .AND. (LB2.GE.-17.AND.LB2.LE.-14)) GOTO 5699
       IF (LB2.EQ.21 .AND. (LB1.GE.-17.AND.LB1.LE.-14)) GOTO 5699
       IF( (((LB1.eq.1.or.LB1.eq.2).or.(LB1.ge.6.and.LB1.le.13))
     &       .AND.(LB2.GE.-17.AND.LB2.LE.-14)) .OR.
     &     (((LB2.eq.1.or.LB2.eq.2).or.(LB2.ge.6.and.LB2.le.13))
     &      .AND.(LB1.GE.-17.AND.LB1.LE.-14)) )go to 5999
       IF( (((LB1.eq.-1.or.LB1.eq.-2).or.(LB1.le.-6.and.LB1.ge.-13))
     &       .AND.(LB2.GE.14.AND.LB2.LE.17)) .OR.
     &     (((LB2.eq.-1.or.LB2.eq.-2).or.(LB2.le.-6.and.LB2.ge.-13))
     &       .AND.(LB1.GE.14.AND.LB1.LE.17)) )go to 5999 
       if(lb1.eq.21.and.lb2.eq.23) go to 8699
       if(lb2.eq.21.and.lb1.eq.23) go to 8699
       if(lb1.eq.30.and.lb2.eq.21) go to 8699
       if(lb2.eq.30.and.lb1.eq.21) go to 8699
       if(lb1.eq.-30.and.lb2.eq.23) go to 8699
       if(lb2.eq.-30.and.lb1.eq.23) go to 8699
       if(lb1.eq.-30.and.lb2.eq.30) go to 8699
       if(lb2.eq.-30.and.lb1.eq.30) go to 8699
       IF( ((lb1.eq.23.or.lb1.eq.21.or.iabs(lb1).eq.30) .and.
     &      (lb2.ge.25.and.lb2.le.28)) .OR.
     &     ((lb2.eq.23.or.lb2.eq.21.or.iabs(lb2).eq.30) .and.
     &      (lb1.ge.25.and.lb1.le.28)) ) go to 8799
       IF( (iabs(lb1).eq.30.and.(lb2.ge.3.and.lb2.le.5)) .OR.
     &     (iabs(lb2).eq.30.and.(lb1.ge.3.and.lb1.le.5)) )go to 8799
       IF( (lb1.eq.29 .and.(lb2.eq.1.or.lb2.eq.2.or.
     &       (lb2.ge.6.and.lb2.le.9))) .OR.
     &     (lb2.eq.29 .and.(lb1.eq.1.or.lb1.eq.2.or.
     &       (lb1.ge.6.and.lb1.le.9))) )go to 7222
       IF( (lb1.eq.29 .and.((lb2.ge.3.and.lb2.le.5).or.
     &      (lb2.ge.21.and.lb2.le.28).or.iabs(lb2).eq.30)) .OR.
     &     (lb2.eq.29 .and.((lb1.ge.3.and.lb1.le.5).or.
     &      (lb1.ge.21.and.lb1.le.28).or.iabs(lb1).eq.30)) )THEN
             go to 7444
      endif
      if( ((iabs(lb1).ge.14.and.iabs(lb1).le.17).or.iabs(lb1).ge.40)
     &    .and.((lb2.ge.25.and.lb2.le.29).or.lb2.eq.0) )go to 888
      if( ((iabs(lb2).ge.14.and.iabs(lb2).le.17).or.iabs(lb2).ge.40)
     &    .and.((lb1.ge.25.and.lb1.le.29).or.lb1.eq.0) )go to 888
        if( ((lb1.eq.23.or.lb1.eq.30).and.(lb2.eq.1.or.lb2.eq.2.or.
     &         (lb2.ge.6.and.lb2.le.13))) .OR.
     &      ((lb2.eq.23.or.lb2.eq.30).and.(lb1.eq.1.or.lb1.eq.2.or.
     &         (lb1.ge.6.and.lb1.le.13))) ) go to 888
        if( ((lb1.eq.21.or.lb1.eq.-30).and.(lb2.eq.-1.or.lb2.eq.-2.or.
     &       (lb2.ge.-13.and.lb2.le.-6))) .OR. 
     &      ((lb2.eq.21.or.lb2.eq.-30).and.(lb1.eq.-1.or.lb1.eq.-2.or.
     &       (lb1.ge.-13.and.lb1.le.-6))) ) go to 888
       If( ((lb1.ge.14.and.lb1.le.17).and.(lb2.ge.6.and.lb2.le.13))
     & .OR.((lb2.ge.14.and.lb2.le.17).and.(lb1.ge.6.and.lb1.le.13)) )
     &   go to 7799
       If(((lb1.le.-14.and.lb1.ge.-17).and.(lb2.le.-6.and.lb2.ge.-13))
     &.OR.((lb2.le.-14.and.lb2.ge.-17).and.(lb1.le.-6.and.lb1.ge.-13)))
     &   go to 7799
       if( iabs(lb1).ge.40 .or. iabs(lb2).ge.40
     &    .or. (lb1.le.-14.and.lb1.ge.-17) 
     &    .or. (lb2.le.-14.and.lb2.ge.-17) )go to 400
        if((lb1.eq.-1.or.lb1.eq.-2.or.(lb1.ge.-13.and.lb1.le.-6))
     1 .and.(lb2.eq.1.or.lb2.eq.2.or.(lb2.ge.6.and.lb2.le.13))) then
            GOTO 2799
       else if((lb2.eq.-1.or.lb2.eq.-2.or.(lb2.ge.-13.and.lb2.le.-6))
     1 .and.(lb1.eq.1.or.lb1.eq.2.or.(lb1.ge.6.and.lb1.le.13))) then
            GOTO 2799
         END IF
         inewka=irun
        call newka(icase,inewka,iseed,dt,nt,
     &                  ictrl,i1,i2,srt,pcx,pcy,pcz,iblock)
        IF (ICTRL .EQ. 1) GOTO 400
       if((iabs(lb1).ge.14.and.iabs(lb1).le.17).
     &  or.(iabs(lb2).ge.14.and.iabs(lb2).le.17))go to 400
       IF((lb1.ge.3.and.lb1.le.5).and.(lb2.ge.3.and.lb2.le.5))GO TO 777
       if(lb1.eq.0.and.(lb2.ge.3.and.lb2.le.5)) go to 777
       if(lb2.eq.0.and.(lb1.ge.3.and.lb1.le.5)) go to 777
       if(lb1.eq.0.and.lb2.eq.0)go to 777
       if( (lb1.ge.25.and.lb1.le.28).and.
     &     (lb2.ge.25.and.lb2.le.28) )goto 777
      If((lb1.ge.25.and.lb1.le.28).and.(lb2.ge.3.and.lb2.le.5))go to 777
      If((lb2.ge.25.and.lb2.le.28).and.(lb1.ge.3.and.lb1.le.5))go to 777
       if((lb1.ge.25.and.lb1.le.28).and.lb2.eq.0)go to 777
       if((lb2.ge.25.and.lb2.le.28).and.lb1.eq.0)go to 777
       if((lb1.eq.23.or.lb1.eq.21).and.(lb2.ge.3.and.lb2.le.5))go to 889
       if((lb2.eq.23.or.lb2.eq.21).and.(lb1.ge.3.and.lb1.le.5))go to 889
        If(iabs(lb1).eq.30.or.iabs(lb2).eq.30) go to 400
        If(lb1.eq.21.or.lb2.eq.21) go to 400
        If(lb1.eq.23.or.lb2.eq.23) go to 400
           IF( (LB1.ge.3.and.LB1.le.5) .and. 
     &         (iabs(LB2).eq.1.or.iabs(LB2).eq.2.or.
     &          (iabs(LB2).ge.6.and.iabs(LB2).le.13)) )GO TO 3
           IF( (LB2.ge.3.and.LB2.le.5) .and. 
     &         (iabs(LB1).eq.1.or.iabs(LB1).eq.2.or.
     &          (iabs(LB1).ge.6.and.iabs(LB1).le.13)) )GO TO 3
           IF( (LB1.ge.25.and.LB1.le.28) .and. 
     &         (iabs(LB2).eq.1.or.iabs(LB2).eq.2.or.
     &          (iabs(LB2).ge.6.and.iabs(LB2).le.13)) )GO TO 33
           IF( (LB2.ge.25.and.LB2.le.28) .and. 
     &         (iabs(LB1).eq.1.or.iabs(LB1).eq.2.or.
     &          (iabs(LB1).ge.6.and.iabs(LB1).le.13)) )GO TO 33
           IF( LB1.eq.0 .and. 
     &         (iabs(LB2).eq.1.or.iabs(LB2).eq.2.or.
     &          (iabs(LB2).ge.6.and.iabs(LB2).le.13)) )GO TO 547
           IF( LB2.eq.0 .and. 
     &         (iabs(LB1).eq.1.or.iabs(LB1).eq.2.or.
     &          (iabs(LB1).ge.6.and.iabs(LB1).le.13)) )GO TO 547
            IF((LB1.eq.1.or.lb1.eq.2).
     &        AND.(LB2.GT.5.and.lb2.le.13))GOTO 44
            IF((LB2.eq.1.or.lb2.eq.2).
     &        AND.(LB1.GT.5.and.lb1.le.13))GOTO 44
            IF((LB1.eq.-1.or.lb1.eq.-2).
     &        AND.(LB2.LT.-5.and.lb2.ge.-13))GOTO 44
            IF((LB2.eq.-1.or.lb2.eq.-2).
     &        AND.(LB1.LT.-5.and.lb1.ge.-13))GOTO 44
       IF((LB1.eq.1.or.lb1.eq.2).AND.(LB2.eq.1.or.lb2.eq.2))GOTO 4
       IF((LB1.eq.-1.or.lb1.eq.-2).AND.(LB2.eq.-1.or.lb2.eq.-2))GOTO 4
            IF((LB1.GT.5.and.lb1.le.13).AND.
     &         (LB2.GT.5.and.lb2.le.13)) GOTO 444
            IF((LB1.LT.-5.and.lb1.ge.-13).AND.
     &         (LB2.LT.-5.and.lb2.ge.-13)) GOTO 444
       if((lb1.lt.3).and.(lb2.ge.14.and.lb2.le.17))goto 400
       if((lb2.lt.3).and.(lb1.ge.14.and.lb1.le.17))goto 400
       if((lb1.ge.14.and.lb1.le.17).and.
     &  (lb2.ge.14.and.lb2.le.17))goto 400
              go to 400
547           IF(LB1*LB2.EQ.0)THEN
           ece=(em1+em2+0.02)**2
           xkaon0=0.
           if(srt.ge.1.63.AND.SRT.LE.1.7)xkaon0=pnlka(srt)
           IF(SRT.GT.1.7)XKAON0=PNLKA(SRT)+pnska(srt)
            XKAON0 = 2.0 * XKAON0
           xkaon=xkaon0
            XETA=XN1535(I1,I2,0)
        If((iabs(LB(I1)).ge.6.and.iabs(LB(I1)).le.13).or.
     &     (iabs(LB(I2)).ge.6.and.iabs(LB(I2)).le.13)) xeta=0.      
            IF((XETA+xkaon).LE.1.e-06)GO TO 400
            DSE=SQRT((XETA+XKAON)/PI)
           DELTRE=DSE+0.1
        px1cm=pcx
        py1cm=pcy
        pz1cm=pcz
         CALL DISTCE(I1,I2,DELTRE,DSE,DT,ECE,SRT,IC,
     1   PCX,PCY,PCZ)
         IF(IC.EQ.-1) GO TO 400
         ekaon(4,iss)=ekaon(4,iss)+1
        IF(XKAON0/(XKAON+XETA).GT.RANART(NSEED))then
        CALL CREN(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK)
       IF(IBLOCK.EQ.7) then
          LPN=LPN+1
       elseIF(IBLOCK.EQ.-7) then
       endif
       em1=e(i1)
       em2=e(i2)
       GO TO 440
       endif
        resona=1.
         GO TO 98
         ENDIF
3           CONTINUE
           px1cm=pcx
           py1cm=pcy
           pz1cm=pcz
           xkaon0=0.
           if(srt.ge.1.63.AND.SRT.LE.1.7)xkaon0=pnlka(srt)
           IF(SRT.GT.1.7)XKAON0=PNLKA(SRT)+pnska(srt)
            XKAON0 = 2.0 * XKAON0
         Xphi = 0.
       if( ( ((lb1.ge.1.and.lb1.le.2).or.
     &        (lb1.ge.6.and.lb1.le.9))
     &   .OR.((lb2.ge.1.and.lb2.le.2).or.
     &        (lb2.ge.6.and.lb2.le.9)) )
     &       .AND. srt.gt.1.958)
     &        call pibphi(srt,lb1,lb2,em1,em2,Xphi,xphin)
         If((iabs(LB(I1)).ge.6.and.iabs(LB(I1)).le.13).or.
     &      (iabs(LB(I2)).ge.6.and.iabs(LB(I2)).le.13)) go to 31
            EC=(em1+em2+0.02)**2
           xkaon=0.
           if(srt.gt.1.23)xkaon=(pionpp(srt)+PIPP1(SRT))/2.
           IF( (LB1*LB2.EQ.5.OR.((LB1*LB2.EQ.6).AND.
     &         (LB1.EQ.3.OR.LB2.EQ.3)))
     &    .OR. (LB1*LB2.EQ.-3.OR.((LB1*LB2.EQ.-10).AND.
     &         (LB1.EQ.5.OR.LB2.EQ.5))) )then    
              XMAX=190.
              xmaxn=0
              xmaxn1=0
              xdirct=dirct1(srt)
               go to 678
           endif
            IF( (LB1*LB2.EQ.3.OR.((LB1*LB2.EQ.10).AND.
     &          (LB1.EQ.5.OR.LB2.EQ.5)))
     &     .OR. (LB1*LB2.EQ.-5.OR.((LB1*LB2.EQ.-6).AND.
     &          (LB1.EQ.3.OR.LB2.EQ.3))) )then      
              XMAX=27.
              xmaxn=2./3.*25.*0.6
               xmaxn1=2./3.*40.*0.5
              xdirct=dirct2(srt)
               go to 678
              endif
            IF((LB1.EQ.4.OR.LB2.EQ.4).AND.
     &         (iabs(LB1*LB2).EQ.4.OR.iabs(LB1*LB2).EQ.8))then
              XMAX=50.
              xmaxn=1./3.*25*0.6
              xmaxn1=1/3.*40.*0.5
              xdirct=dirct3(srt)
                go to 678
              endif
678           xnpin1=0
           xnpin=0
            XNPID=XNPI(I1,I2,1,XMAX)
           if(xmaxn1.ne.0)xnpin1=XNPI(i1,i2,2,XMAXN1)
            if(xmaxn.ne.0)XNPIN=XNPI(I1,I2,0,XMAXN)
           xres=xnpid+xnpin+xnpin1
           xnelas=xres+xdirct 
           icheck=1
           go to 34
31           ec=(em1+em2+0.02)**2
           xreab=reab(i1,i2,srt,1)
          if((iabs(lb1).ge.10.and.iabs(lb1).le.13)
     1         .or.(iabs(lb2).ge.10.and.iabs(lb2).le.13)) XREAB = 0.
           xkaon=xkaon0+xreab
        IF((iabs(LB1).GT.9.AND.iabs(LB1).LE.13) .OR.
     &      (iabs(LB2).GT.9.AND.iabs(LB2).LE.13))THEN
           Xnelas=1.0
        ELSE
           XNELAS=DPION(EM1,EM2,LB1,LB2,SRT)
        ENDIF
           icheck=2
34          IF((Xnelas+xkaon+Xphi).LE.0.000001)GO TO 400
            DS=SQRT((Xnelas+xkaon+Xphi)/PI)
           deltar=ds+0.1
         CALL DISTCE(I1,I2,DELTAR,DS,DT,EC,SRT,IC,
     1   PCX,PCY,PCZ)
         IF(IC.EQ.-1) GO TO 400
       ekaon(4,iss)=ekaon(4,iss)+1
        if(icheck.eq.2)then
      if(xnelas/(xnelas+xkaon+Xphi).ge.RANART(NSEED))then
               call Crdir(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK)
              go to 440
              else
               go to 96
                endif
              endif
       IF((XKAON+Xphi)/(XKAON+Xphi+Xnelas).GT.RANART(NSEED))GO TO 95
        if(xdirct/xnelas.ge.RANART(NSEED))then
               call Crdir(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK)
              go to 440
              endif
           IF( (LB1*LB2.EQ.5.OR.((LB1*LB2.EQ.6).AND.
     &         (LB1.EQ.3.OR.LB2.EQ.3)))
     &    .OR. (LB1*LB2.EQ.-3.OR.((LB1*LB2.EQ.-10).AND.
     &         (LB1.EQ.5.OR.LB2.EQ.5))) )then    
        GO TO 99
       else
            XX=(XNPIN+xnpin1)/xres
            IF(RANART(NSEED).LT.XX)THEN
        xx0=xnpin/(xnpin+xnpin1)
        if(RANART(NSEED).lt.xx0)then
         RESONA=0.
         GO TO 97
        else
        resona=1.
         GO TO 98
        endif
         ELSE
         GO TO 99
         ENDIF
         ENDIF
97       CONTINUE
            IF(RESONA.EQ.0.)THEN
            I=I1
            IF(EM1.LT.0.6)I=I2
           IF( (LB1*LB2.EQ.10.AND.(LB1.EQ.5.OR.LB2.EQ.5))
     &      .OR.(LB1*LB2.EQ.-6.AND.(LB1.EQ.3.OR.LB2.EQ.3)) )THEN
            LB(I)=11
           go to 303
            ENDIF
            IF(iabs(LB(I1)*LB(I2)).EQ.4.AND.
     &         (LB(I1).EQ.4.OR.LB(I2).EQ.4))THEN    
            LB(I)=11
           go to 303
            ENDIF
            IF(iabs(LB(I1)*LB(I2)).EQ.8.AND.
     &        (LB(I1).EQ.4.OR.LB(I2).EQ.4))THEN    
            LB(I)=10
           go to 303
            ENDIF
            IF( (LB(I1)*LB(I2).EQ.3)
     &      .OR.(LB(I1)*LB(I2).EQ.-5) )THEN
            LB(I)=10
            ENDIF
303         CALL DRESON(I1,I2)
            if(LB1.lt.0.or.LB2.lt.0) LB(I)=-LB(I)
            lres=lres+1
            GO TO 101
            ENDIF
98          IF(RESONA.EQ.1.)THEN
            I=I1
            IF(EM1.LT.0.6)I=I2
            IF( (LB1*LB2.EQ.10.AND.(LB1.EQ.5.OR.LB2.EQ.5))
     &      .OR.(LB1*LB2.EQ.-6.AND.(LB1.EQ.3.OR.LB2.EQ.3)) )THEN
            LB(I)=13
           go to 304
            ENDIF
            IF(iabs(LB(I1)*LB(I2)).EQ.4.AND.
     &           (LB(I1).EQ.4.OR.LB(I2).EQ.4))THEN 
            LB(I)=13
           go to 304
            ENDIF
            IF(iabs(LB(I1)*LB(I2)).EQ.8.AND.
     &           (LB(I1).EQ.4.OR.LB(I2).EQ.4))THEN      
            LB(I)=12
           go to 304
            ENDIF
            IF( (LB(I1)*LB(I2).EQ.3)
     &      .OR.(LB(I1)*LB(I2).EQ.-5) )THEN
            LB(I)=12
           go to 304
           endif
           if(lb(i1)*lb(i2).eq.0)then
            if(iabs(lb(i1)).eq.1.or.iabs(lb(i2)).eq.1)then
           LB(I)=13
           go to 304
           ELSE
           LB(I)=12
           ENDIF
           endif
304         CALL DRESON(I1,I2)
            if(LB1.lt.0.or.LB2.lt.0) LB(I)=-LB(I) 
            lres=lres+1
            GO TO 101
            ENDIF
99      LRES=LRES+1
        I=I1
        IF(EM1.LE.0.6)I=I2
            IF( (LB(I1)*LB(I2).EQ.5)
     &      .OR.(LB(I1)*LB(I2).EQ.-3) )THEN
        LB(I)=9
       go to 305
        ENDIF
       IF(iabs(LB(I1)*LB(I2)).EQ.4.AND.(LB(I1).EQ.4.OR.LB(I2).EQ.4))then
        LB(I)=8
       go to 305
        ENDIF
       IF( (LB(I1)*LB(I2).EQ.10.AND.(LB(I1).EQ.5.OR.LB(I2).EQ.5))
     & .OR.(LB(I1)*LB(I2).EQ.-6.AND.(LB(I1).EQ.3.OR.LB(I2).EQ.3)) )THEN
        LB(I)=8
       go to 305
        ENDIF
       IF(iabs(LB(I1)*LB(I2)).EQ.8.AND.(LB(I1).EQ.4.OR.LB(I2).EQ.4))THEN
        LB(I)=7
       go to 305
        ENDIF
            IF( (LB(I1)*LB(I2).EQ.3)
     &      .OR.(LB(I1)*LB(I2).EQ.-5) )THEN
        LB(I)=7
       go to 305
        ENDIF
       IF( (LB(I1)*LB(I2).EQ.6.AND.(LB(I1).EQ.3.OR.LB(I2).EQ.3))
     & .OR.(LB(I1)*LB(I2).EQ.-10.AND.(LB(I1).EQ.5.OR.LB(I2).EQ.5)) )THEN 
        LB(I)=6
        ENDIF
305     CALL DRESON(I1,I2)
        if(LB1.lt.0.or.LB2.lt.0) LB(I)=-LB(I) 
       GO TO 101
889       CONTINUE
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=(em1+em2+0.02)**2
       spika=60./(1.+4.*(srt-0.895)**2/(0.05)**2)
        call Crkpla(PX1CM,PY1CM,PZ1CM,EC,SRT,spika,
     &                  emm1,emm2,lbp1,lbp2,I1,I2,icase,srhoks)
         if(icase .eq. 0) then
            iblock=0
            go to 400
         endif
       if(icase .eq. 1)then
             call KSRESO(I1,I2)
             iblock = 171
              lres=lres+1
              go to 101
       elseif(icase .eq. 2)then
             iblock = 71
       elseif(iabs(icase).eq.5)then
             iblock = 88
       else
             iblock = 222
       endif
             LB(I1) = lbp1
             LB(I2) = lbp2
             E(I1) = emm1
             E(I2) = emm2
             em1=e(i1)
             em2=e(i2)
             ntag = 0
             go to 440
33       continue
       em1=e(i1)
       em2=e(i2)
           xelstc=0
            if((lb1.ge.25.and.lb1.le.28).and.
     &    (iabs(lb2).eq.1.or.iabs(lb2).eq.2))
     &      xelstc=ERHON(EM1,EM2,LB1,LB2,SRT)
            if((lb2.ge.25.and.lb2.le.28).and.
     &   (iabs(lb1).eq.1.or.iabs(lb1).eq.2))
     &      xelstc=ERHON(EM1,EM2,LB1,LB2,SRT)
            ec=(em1+em2+0.02)**2
           xkaon0=0
           if(srt.ge.1.63.AND.SRT.LE.1.7)xkaon0=pnlka(srt)
           IF(SRT.GT.1.7)XKAON0=PNLKA(SRT)+pnska(srt)
           if(xkaon0.lt.0)xkaon0=0
            XKAON0 = 2.0 * XKAON0
           xkaon=xkaon0
           ichann=0
         Xphi = 0.
       if( ( (((lb1.ge.1.and.lb1.le.2).or.
     &         (lb1.ge.6.and.lb1.le.9))
     &         .and.(lb2.ge.25.and.lb2.le.27))
     &   .OR.(((lb2.ge.1.and.lb2.le.2).or.
     &         (lb2.ge.6.and.lb2.le.9))
     &        .and.(lb1.ge.25.and.lb1.le.27)) ).AND. srt.gt.1.958)
     &    call pibphi(srt,lb1,lb2,em1,em2,Xphi,xphin)
        if((iabs(lb1).ge.6.and.lb2.ge.25).or.
     &    (lb1.ge.25.and.iabs(lb2).ge.6))then
           ichann=1
           ictrl=2
           if(lb1.eq.28.or.lb2.eq.28)ictrl=3
            xreab=reab(i1,i2,srt,ictrl)
            if((iabs(lb1).ge.10.and.iabs(lb1).le.13)
     1           .or.(iabs(lb2).ge.10.and.iabs(lb2).le.13)) XREAB = 0.
        if(xreab.lt.0)xreab=1.E-06
            xkaon=xkaon0+xreab
          XELSTC=1.0
           endif
            DS=SQRT((XKAON+Xphi+xelstc)/PI)
        DELTAR=DS+0.1
       px1cm=pcx
       py1cm=pcy
       pz1cm=pcz
         CALL DISTCE(I1,I2,DELTAR,DS,DT,EC,SRT,IC,
     1   PCX,PCY,PCZ)
         IF(IC.EQ.-1) GO TO 400
        ekaon(4,iss)=ekaon(4,iss)+1
       if(xelstc/(xelstc+xkaon+Xphi).gt.RANART(NSEED))then
       call crdir(px1CM,py1CM,pz1CM,srt,I1,i2,IBLOCK)
       go to 440
       endif
        CALL CRRD(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     1  IBLOCK,xkaon0,xkaon,Xphi,xphin)
       IF(IBLOCK.EQ.7) then
          LPN=LPN+1
       elseIF(IBLOCK.EQ.-7) then
       endif
       if(iblock.eq.81) lrhor=lrhor+1
       if(iblock.eq.82) lomgar=lomgar+1
       em1=e(i1)
       em2=e(i2)
       GO TO 440
95       continue
        CALL CRPN(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     1  IBLOCK,xkaon0,xkaon,Xphi,xphin)
       IF(IBLOCK.EQ.7) then
          LPN=LPN+1
       elseIF(IBLOCK.EQ.-7) then
       endif
       if(iblock.eq.77) lpd=lpd+1
       if(iblock.eq.78) lrho=lrho+1
       if(iblock.eq.79) lomega=lomega+1
       em1=e(i1)
       em2=e(i2)
       GO TO 440
96       continue
        CALL CRPD(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     1  IBLOCK,xkaon0,xkaon,Xphi,xphin)
       IF(IBLOCK.EQ.7) then
          LPN=LPN+1
       elseIF(IBLOCK.EQ.-7) then
       endif
       if(iblock.eq.80) lpdr=lpdr+1
       em1=e(i1)
       em2=e(i2)
       GO TO 440
101       continue
        IF(E(I2).EQ.0.)GO TO 600
        IF(E(I1).EQ.0.)GO TO 800
44      CONTINUE
       cutoff=em1+em2+0.02
       IF(SRT.LE.CUTOFF)GO TO 400
        IF(SRT.GT.2.245)THEN
       SIGNN=PP2(SRT)
       ELSE
        SIGNN = 35.0 / (1. + (SRT - CUTOFF) * 100.0)  +  20.0
       ENDIF 
        call XND(pcx,pcy,pcz,srt,I1,I2,xinel,
     &               sigk,xsk1,xsk2,xsk3,xsk4,xsk5)
       sig=signn+xinel
        EC=(EM1+EM2+0.02)**2
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        ianti=0
        if(lb(i1).lt.0 .and. lb(i2).lt.0) ianti=1
        call sbbdm(srt,sdprod,ianti,lbm,xmm,pfinal)
        sig=sig+sdprod
        ipdflag=0
        if(idpert.eq.1) then
           ipert1=1
           sigr0=sig
           dspert=sqrt(sigr0/pi/10.)
           dsrpert=dspert+0.1
           CALL DISTCE(I1,I2,dsrpert,dspert,DT,EC,SRT,IC,
     1          PX1CM,PY1CM,PZ1CM)
           IF(IC.EQ.-1) GO TO 363
           signn0=0.
           CALL CRND(IRUN,PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     &  IBLOCK,SIGNN0,SIGr0,sigk,xsk1,xsk2,xsk3,xsk4,xsk5,NT,ipert1)
           ipdflag=1
 363       continue
           ipert1=0
        endif
        if(idpert.eq.2) ipert1=1
        DS=SQRT(SIG/(10.*PI))
        DELTAR=DS+0.1
        CALL DISTCE(I1,I2,DELTAR,DS,DT,EC,SRT,IC,
     1  PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) then
           if(ipdflag.eq.1) iblock=501
           GO TO 400
        endif
        ekaon(3,iss)=ekaon(3,iss)+1
        go to 361
 361    continue 
        CALL CRND(IRUN,PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     &     IBLOCK,SIGNN,SIG,sigk,xsk1,xsk2,xsk3,xsk4,xsk5,NT,ipert1)
        IF(iblock.eq.0.and.ipdflag.eq.1) iblock=501
        IF(IBLOCK.EQ.11)THEN
           LNDK=LNDK+1
           GO TO 400
        elseIF(IBLOCK.EQ.-11.or.iblock.eq.501) then
           GO TO 400
        ENDIF
        if(iblock .eq. 222)then
           GO TO 400
        ENDIF
        em1=e(i1)
        em2=e(i2)
        GO TO 440
4       CONTINUE
        CUTOFF=em1+em2+0.14
        IF(SRT.GT.2.245)THEN
           SIG=ppt(srt)
           SIGNN=SIG-PP1(SRT)
        ELSE
           SIG=XPP(SRT)
           IF(ZET(LB(I1))*ZET(LB(I2)).LE.0)SIG=XNP(SRT)
           IF(ZET(LB(I1))*ZET(LB(I2)).GT.0)SIG=XPP(SRT)
           IF(ZET(LB(I1)).EQ.0.
     &          AND.ZET(LB(I2)).EQ.0)SIG=XPP(SRT)
           if((lb(i1).eq.-1.and.lb(i2).eq.-2) .or.
     &          (lb(i2).eq.-1.and.lb(i1).eq.-2))sig=xnp(srt)
           IF (SRT .LT. 1.897) THEN
              SIGNN = SIG
           ELSE 
              SIGNN = 35.0 / (1. + (SRT - 1.897) * 100.0)  +  20.0
           ENDIF
        ENDIF 
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        ianti=0
        if(lb(i1).lt.0 .and. lb(i2).lt.0) ianti=1
        call sbbdm(srt,sdprod,ianti,lbm,xmm,pfinal)
        sig=sig+sdprod
        ipdflag=0
        if(idpert.eq.1) then
           ipert1=1
           EC=2.012**2
           sigr0=sig
           dspert=sqrt(sigr0/pi/10.)
           dsrpert=dspert+0.1
           CALL DISTCE(I1,I2,dsrpert,dspert,DT,EC,SRT,IC,
     1          PX1CM,PY1CM,PZ1CM)
           IF(IC.EQ.-1) GO TO 365
           signn0=0.
           CALL CRNN(IRUN,PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK,
     1          NTAG,signn0,sigr0,NT,ipert1)
           ipdflag=1
 365       continue
           ipert1=0
        endif
        if(idpert.eq.2) ipert1=1
        IF(SIGNN.LE.0) then
           if(ipdflag.eq.1) iblock=501
           GO TO 400
        endif
        EC=3.59709
        ds=sqrt(sig/pi/10.)
        dsr=ds+0.1
        IF((E(I1).GE.1.).AND.(e(I2).GE.1.))EC=4.75
        CALL DISTCE(I1,I2,dsr,ds,DT,EC,SRT,IC,
     1       PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) then
           if(ipdflag.eq.1) iblock=501
           GO TO 400
        endif
        go to 362
 362    ekaon(1,iss)=ekaon(1,iss)+1
        CALL CRNN(IRUN,PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK,
     1       NTAG,SIGNN,SIG,NT,ipert1)
        IF(iblock.eq.0.and.ipdflag.eq.1) iblock=501
        IF(IBLOCK.EQ.4.OR.IBLOCK.Eq.9.or.iblock.ge.44.OR.IBLOCK.EQ.-9
     &       .or.iblock.eq.222.or.iblock.eq.501)THEN
           LCOLL=LCOLL+1
           if(iblock.eq.4)then
              LDIRT=LDIRT+1
           elseif(iblock.eq.44)then
              LDdrho=LDdrho+1
           elseif(iblock.eq.45)then
              Lnnrho=Lnnrho+1
           elseif(iblock.eq.46)then
              Lnnom=Lnnom+1
           elseif(iblock .eq. 222)then
           elseIF(IBLOCK.EQ.9) then
              LNNK=LNNK+1
           elseIF(IBLOCK.EQ.-9) then
           endif
           GO TO 400
        ENDIF
        em1=e(i1)
        em2=e(i2)
        GO TO 440
 505    continue
        ianti=0
        if(lb(i1).lt.0 .or. lb(i2).lt.0) ianti=1
        call sdmbb(SRT,sdm,ianti)
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=2.012**2
        ds=sqrt(sdm/31.4)
        dsr=ds+0.1
        CALL DISTCE(I1,I2,dsr,ds,DT,EC,SRT,IC,PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
        CALL crdmbb(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK,
     1       NTAG,sdm,NT,ianti)
        LCOLL=LCOLL+1
        GO TO 400
 506    continue
        ianti=0
        if(lb(i1).lt.0 .or. lb(i2).lt.0) ianti=1
        call sdbelastic(SRT,sdb)
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=2.012**2
        ds=sqrt(sdb/31.4)
        dsr=ds+0.1
        CALL DISTCE(I1,I2,dsr,ds,DT,EC,SRT,IC,PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
        CALL crdbel(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK,
     1       NTAG,sdb,NT,ianti)
        LCOLL=LCOLL+1
        GO TO 400
 444    CONTINUE
       CUTOFF=em1+em2+0.02
       IF(SRT.LE.CUTOFF)GO TO 400
        IF(SRT.GT.2.245)THEN
       SIGNN=PP2(SRT)
       ELSE
        SIGNN = 35.0 / (1. + (SRT - CUTOFF) * 100.0)  +  20.0
       ENDIF 
       IF(SIGNN.LE.0)GO TO 400
      CALL XDDIN(PCX,PCY,PCZ,SRT,I1,I2,
     &XINEL,SIGK,XSK1,XSK2,XSK3,XSK4,XSK5)
       SIG=SIGNN+XINEL
       EC=(EM1+EM2+0.02)**2
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        ianti=0
        if(lb(i1).lt.0 .and. lb(i2).lt.0) ianti=1
        call sbbdm(srt,sdprod,ianti,lbm,xmm,pfinal)
        sig=sig+sdprod
        ipdflag=0
        if(idpert.eq.1) then
           ipert1=1
           sigr0=sig
           dspert=sqrt(sigr0/pi/10.)
           dsrpert=dspert+0.1
           CALL DISTCE(I1,I2,dsrpert,dspert,DT,EC,SRT,IC,
     1          PX1CM,PY1CM,PZ1CM)
           IF(IC.EQ.-1) GO TO 367
           signn0=0.
           CALL CRDD(IRUN,PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     1          IBLOCK,NTAG,SIGNN0,SIGr0,NT,ipert1)
           ipdflag=1
 367       continue
           ipert1=0
        endif
        if(idpert.eq.2) ipert1=1
        ds=sqrt(sig/31.4)
        dsr=ds+0.1
        CALL DISTCE(I1,I2,dsr,ds,DT,EC,SRT,IC,
     1  PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) then
           if(ipdflag.eq.1) iblock=501
           GO TO 400
        endif
       go to 364
364       ekaon(2,iss)=ekaon(2,iss)+1
        CALL CRDD(IRUN,PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     1  IBLOCK,NTAG,SIGNN,SIG,NT,ipert1)
        IF(iblock.eq.0.and.ipdflag.eq.1) iblock=501
        IF(iabs(IBLOCK).EQ.10)THEN
           LCOLL=LCOLL+1
           IF(IBLOCK.EQ.10)THEN
              LDDK=LDDK+1
           elseIF(IBLOCK.EQ.-10) then
           endif
           GO TO 400
        ENDIF
        if(iblock .eq. 222.or.iblock.eq.501)then
           GO TO 400
        ENDIF
        em1=e(i1)
        em2=e(i2)
        GO TO 440
777       CONTINUE
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
       ec0=em1+em2+0.02
       IF(SRT.LE.ec0)GO TO 400
       ec=(em1+em2+0.02)**2
       ppel=20.
        ipp=1
       if(lb1.lt.3.or.lb1.gt.5.or.lb2.lt.3.or.lb2.gt.5)go to 778       
       CALL PPXS(LB1,LB2,SRT,PPSIG,spprho,IPP)
       ppel=ppsig
778       ppink=pipik(srt)
        ppink = 2.0 * ppink
       if(lb1.ge.25.and.lb2.ge.25) ppink=rrkk
        if( ( (lb1.eq.0.or.(lb1.ge.3.and.lb1.le.5))
     1       .and.(lb2.ge.25.and.lb2.le.28))
     2       .or. ( (lb2.eq.0.or.(lb2.ge.3.and.lb2.le.5))
     3       .and.(lb1.ge.25.and.lb1.le.28))) then
           ppink=0.
           if(srt.ge.(aka+aks)) ppink = prkk
        endif
        call spprr(lb1,lb2,srt)
        call sppee(lb1,lb2,srt)
        call spppe(lb1,lb2,srt)
        call srpre(lb1,lb2,srt)
        call sopoe(lb1,lb2,srt)
        call srree(lb1,lb2,srt)
        ppinnb=0.
        if(srt.gt.thresh(1)) then
           call getnst(srt)
           if(lb1.ge.3.and.lb1.le.5.and.lb2.ge.3.and.lb2.le.5) then
              ppinnb=ppbbar(srt)
           elseif((lb1.ge.3.and.lb1.le.5.and.lb2.ge.25.and.lb2.le.27)
     1 .or.(lb2.ge.3.and.lb2.le.5.and.lb1.ge.25.and.lb1.le.27)) then
              ppinnb=prbbar(srt)
           elseif(lb1.ge.25.and.lb1.le.27
     1             .and.lb2.ge.25.and.lb2.le.27) then
              ppinnb=rrbbar(srt)
           elseif((lb1.ge.3.and.lb1.le.5.and.lb2.eq.28)
     1             .or.(lb2.ge.3.and.lb2.le.5.and.lb1.eq.28)) then
              ppinnb=pobbar(srt)
           elseif((lb1.ge.25.and.lb1.le.27.and.lb2.eq.28)
     1             .or.(lb2.ge.25.and.lb2.le.27.and.lb1.eq.28)) then
              ppinnb=robbar(srt)
           elseif(lb1.eq.28.and.lb2.eq.28) then
              ppinnb=oobbar(srt)
           else
              if(lb1.ne.0.and.lb2.ne.0) 
     1             write(6,*) 'missed MM lb1,lb2=',lb1,lb2
           endif
        endif
        ppin=ppink+ppinnb+pprr+ppee+pppe+rpre+xopoe+rree
       if((ppel+ppin).le.0.01)go to 400
       DSPP=SQRT((ppel+ppin)/31.4)
       dsppr=dspp+0.1
        CALL DISTCE(I1,I2,dsppr,DSPP,DT,EC,SRT,IC,
     1  PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
       if(ppel.eq.0)go to 400
       ekaon(5,iss)=ekaon(5,iss)+1
        CALL CRPP(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     1  IBLOCK,ppel,ppin,spprho,ipp)
       if(iblock.eq.666)go to 555
       if(iblock.eq.6)LPP=LPP+1
       if(iblock.eq.66)then
          LPPk=LPPk+1
       elseif(iblock.eq.366)then
          LPPk=LPPk+1
       elseif(iblock.eq.367)then
          LPPk=LPPk+1
       endif
       em1=e(i1)
       em2=e(i2)
       go to 440
2799        CONTINUE
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=(em1+em2+0.02)**2
       DSppb=SQRT(xppbar(srt)/PI/10.)
       dsppbr=dsppb+0.1
        CALL DISTCE(I1,I2,dsppbr,DSppb,DT,EC,SRT,IC,
     1  PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
        CALL Crppba(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     1  IBLOCK)
       em1=e(i1)
       em2=e(i2)
       go to 440
3555    PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=(em1+em2+0.02)**2
       DSkk=SQRT(SIG/PI/10.)
       dskk0=dskk+0.1
        CALL DISTCE(I1,I2,dskk0,DSkk,DT,EC,SRT,IC,
     1  PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
        CALL Crlaba(PX1CM,PY1CM,PZ1CM,SRT,brel,brsgm,
     &                  I1,I2,nt,IBLOCK,nchrg,icase)
       em1=e(i1)
       em2=e(i2)
       go to 440
3455    PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        call pertur(PX1CM,PY1CM,PZ1CM,SRT,IRUN,I1,I2,nt,kp,icontp)
        if(icontp .eq. 0)then
         em1 = e(i1)
         em2 = e(i2)
         iblock = 727
          go to 440
        endif
        if (e(i1) .eq. 0.) go to 800
        if (e(i2) .eq. 0.) go to 600
        go to 400
7222        CONTINUE
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=(em1+em2+0.02)**2
        CALL XphiB(LB1, LB2, EM1, EM2, SRT,
     &             XSK1, XSK2, XSK3, XSK4, XSK5, SIGP)
       DSkk=SQRT(SIGP/PI/10.)
       dskk0=dskk+0.1
        CALL DISTCE(I1,I2,dskk0,DSkk,DT,EC,SRT,IC,
     1  PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
        CALL CRPHIB(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     &     XSK1, XSK2, XSK3, XSK4, XSK5, SIGP, IBLOCK)
       em1=e(i1)
       em2=e(i2)
       go to 440
7444        CONTINUE
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=(em1+em2+0.02)**2
        CALL PHIMES(I1, I2, SRT, XSK1, XSK2, XSK3, XSK4, XSK5,
     1     XSK6, XSK7, SIGPHI)
       DSkk=SQRT(SIGPHI/PI/10.)
       dskk0=dskk+0.1
        CALL DISTCE(I1,I2,dskk0,DSkk,DT,EC,SRT,IC,
     1  PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
        PZRT = p(3,i1)+p(3,i2)
        ER1 = sqrt( p(1,i1)**2+p(2,i1)**2+p(3,i1)**2+E(i1)**2 )
        ER2 = sqrt( p(1,i2)**2+p(2,i2)**2+p(3,i2)**2+E(i2)**2 )
        ERT = ER1+ER2
        yy = 0.5*log( (ERT+PZRT)/(ERT-PZRT) )
        CALL CRPHIM(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     &  XSK1, XSK2, XSK3, XSK4, XSK5, XSK6, SIGPHI, IKKG, IKKL, IBLOCK)
       em1=e(i1)
       em2=e(i2)
       go to 440
 7799    CONTINUE
         PX1CM=PCX
         PY1CM=PCY
         PZ1CM=PCZ
         EC=(em1+em2+0.02)**2
         call lambar(i1,i2,srt,siglab)
        DShn=SQRT(siglab/PI/10.)
        dshnr=dshn+0.1
         CALL DISTCE(I1,I2,dshnr,DShn,DT,EC,SRT,IC,
     1    PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
         CALL Crhb(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK)
        em1=e(i1)
        em2=e(i2)
        go to 440
5699        CONTINUE
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=(em1+em2+0.02)**2
        CALL XKHYPE(I1, I2, SRT, XKY1, XKY2, XKY3, XKY4, XKY5,
     &     XKY6, XKY7, XKY8, XKY9, XKY10, XKY11, XKY12, XKY13,
     &     XKY14, XKY15, XKY16, XKY17, SIGK)
       DSkk=SQRT(sigk/PI)
       dskk0=dskk+0.1
        CALL DISTCE(I1,I2,dskk0,DSkk,DT,EC,SRT,IC,
     1  PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
       if(lb(i1).eq.23 .or. lb(i2).eq.23)then
             IKMP = 1
        else
             IKMP = -1
        endif
        CALL Crkhyp(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     &     XKY1, XKY2, XKY3, XKY4, XKY5,
     &     XKY6, XKY7, XKY8, XKY9, XKY10, XKY11, XKY12, XKY13,
     &     XKY14, XKY15, XKY16, XKY17, SIGK, IKMP,
     1  IBLOCK)
       em1=e(i1)
       em2=e(i2)
       go to 440
5999     CONTINUE
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=(em1+em2+0.02)**2
        sigkp = 15.
        DSkk=SQRT(SIGKP/PI/10.)
        dskk0=dskk+0.1
        CALL DISTCE(I1,I2,dskk0,DSkk,DT,EC,SRT,IC,
     1  PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
        CALL CRLAN(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK)
        em1=e(i1)
        em2=e(i2)
        go to 440
8699     CONTINUE
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=(em1+em2+0.02)**2
         CALL Crkphi(PX1CM,PY1CM,PZ1CM,EC,SRT,IBLOCK,
     &                  emm1,emm2,lbp1,lbp2,I1,I2,ikk,icase,rrkk,prkk)
         if(icase .eq. 0) then
            iblock=0
            go to 400
         endif
         if(lbp1.eq.29.or.lbp2.eq.29) then
        PZRT = p(3,i1)+p(3,i2)
        ER1 = sqrt( p(1,i1)**2+p(2,i1)**2+p(3,i1)**2+E(i1)**2 )
        ER2 = sqrt( p(1,i2)**2+p(2,i2)**2+p(3,i2)**2+E(i2)**2 )
        ERT = ER1+ER2
        yy = 0.5*log( (ERT+PZRT)/(ERT-PZRT) )
             iblock = 222
             ntag = 0
          endif
             LB(I1) = lbp1
             LB(I2) = lbp2
             E(I1) = emm1
             E(I2) = emm2
             em1=e(i1)
             em2=e(i2)
             go to 440
8799     CONTINUE
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=(em1+em2+0.02)**2
         CALL Crksph(PX1CM,PY1CM,PZ1CM,EC,SRT,
     &       emm1,emm2,lbp1,lbp2,I1,I2,ikkg,ikkl,iblock,icase,srhoks)
         if(icase .eq. 0) then
            iblock=0
            go to 400
         endif
         if(lbp1.eq.29.or.lbp2.eq.20) then
        PZRT = p(3,i1)+p(3,i2)
        ER1 = sqrt( p(1,i1)**2+p(2,i1)**2+p(3,i1)**2+E(i1)**2 )
        ER2 = sqrt( p(1,i2)**2+p(2,i2)**2+p(3,i2)**2+E(i2)**2 )
        ERT = ER1+ER2
        yy = 0.5*log( (ERT+PZRT)/(ERT-PZRT) )
          endif
             LB(I1) = lbp1
             LB(I2) = lbp2
             E(I1) = emm1
             E(I2) = emm2
             em1=e(i1)
             em2=e(i2)
             go to 440
888       CONTINUE
        PX1CM=PCX
        PY1CM=PCY
        PZ1CM=PCZ
        EC=(em1+em2+0.02)**2
         sig = 10.
         if(iabs(lb1).eq.14.or.iabs(lb2).eq.14 .or.
     &      iabs(lb1).eq.30.or.iabs(lb2).eq.30)sig=20.
         if(lb1.eq.29.or.lb2.eq.29)sig=5.0
       DSkn=SQRT(sig/PI/10.)
       dsknr=dskn+0.1
        CALL DISTCE(I1,I2,dsknr,DSkn,DT,EC,SRT,IC,
     1  PX1CM,PY1CM,PZ1CM)
        IF(IC.EQ.-1) GO TO 400
        CALL Crkn(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,
     1  IBLOCK)
       em1=e(i1)
       em2=e(i2)
       go to 440
 440    CONTINUE
                 IF(IBLOCK.EQ.0)        GOTO 400
              LCOLL = LCOLL +1
              NTAG = 0
              E1CM    = SQRT (EM1**2 + PX1CM**2 + PY1CM**2 + PZ1CM**2)
              P1BETA  = PX1CM*BETAX + PY1CM*BETAY + PZ1CM*BETAZ
              TRANSF  = GAMMA * ( GAMMA * P1BETA / (GAMMA + 1) + E1CM )
              Pt1I1 = BETAX * TRANSF + PX1CM
              Pt2I1 = BETAY * TRANSF + PY1CM
              Pt3I1 = BETAZ * TRANSF + PZ1CM
              go to 90002
90002              continue
                E2CM    = SQRT (EM2**2 + PX1CM**2 + PY1CM**2 + PZ1CM**2)
                TRANSF  = GAMMA * (-GAMMA*P1BETA / (GAMMA + 1.) + E2CM)
                Pt1I2 = BETAX * TRANSF - PX1CM
                Pt2I2 = BETAY * TRANSF - PY1CM
                Pt3I2 = BETAZ * TRANSF - PZ1CM
              go to 90003
90003           IF(IBLOCK.EQ.1) LCNNE=LCNNE+1
              IF(IBLOCK.EQ.5) LDD=LDD+1
                if(iblock.eq.2) LCNND=LCNND+1
              IF(IBLOCK.EQ.8) LKN=LKN+1
                   if(iblock.eq.43) Ldou=Ldou+1
                IF(IBLOCK.EQ.3) LCNDN=LCNDN+1
              p(1,i1)=pt1i1
              p(2,i1)=pt2i1
              p(3,i1)=pt3i1
              p(1,i2)=pt1i2
              p(2,i2)=pt2i2
              p(3,i2)=pt3i2
                PX1     = P(1,I1)
                PY1     = P(2,I1)
                PZ1     = P(3,I1)
                EM1     = E(I1)
                EM2     = E(I2)
                LB1     = LB(I1)
                LB2     = LB(I2)
                ID(I1)  = 2
                ID(I2)  = 2
                E1      = SQRT( EM1**2 + PX1**2 + PY1**2 + PZ1**2 )
                ID1     = ID(I1)
              go to 90004
90004              continue
            AM1=EM1
            AM2=EM2
  400       CONTINUE
 555        continue
  600     CONTINUE
 798      if(nt.eq.ntmax.and.ipi0dcy.eq.1
     1         .and.i1.eq.(MASSR(IRUN)+MSUM)) then
             do ipion=1,NNN
                if(LPION(ipion,IRUN).eq.4) then
                   wid=7.85e-9
                   call resdec(i1,nt,nnn,wid,idecay,ipion)
                endif
             enddo
          endif
  800   CONTINUE
        N0=MASS+MSUM
        DO 1005 N=N0+1,MASSR(IRUN)+MSUM
        IF(E(N) .GT. 0. .OR. LB(N) .GT. 5000)THEN
        NNN=NNN+1
        RPION(1,NNN,IRUN)=R(1,N)
        RPION(2,NNN,IRUN)=R(2,N)
        RPION(3,NNN,IRUN)=R(3,N)
        if(nt.eq.ntmax) then
           ftpisv(NNN,IRUN)=ftsv(N)
           tfdpi(NNN,IRUN)=tfdcy(N)
        endif
        PPION(1,NNN,IRUN)=P(1,N)
        PPION(2,NNN,IRUN)=P(2,N)
        PPION(3,NNN,IRUN)=P(3,N)
        EPION(NNN,IRUN)=E(N)
        LPION(NNN,IRUN)=LB(N)
        PROPI(NNN,IRUN)=PROPER(N)
        dppion(NNN,IRUN)=dpertp(N)
        ENDIF
 1005 CONTINUE
        MASSRN(IRUN)=NNN+MASS
1000   CONTINUE
        IA=0
        IB=0
        DO 10001 IRUN=1,NUM
        IA=IA+MASSR(IRUN-1)
        IB=IB+MASSRN(IRUN-1)
        DO 10001 IC=1,MASSRN(IRUN)
        IE=IA+IC
        IG=IB+IC
        IF(IC.LE.MASS)THEN
        RT(1,IG)=R(1,IE)
        RT(2,IG)=R(2,IE)
        RT(3,IG)=R(3,IE)
        if(nt.eq.ntmax) then
           fttemp(IG)=ftsv(IE)
           tft(IG)=tfdcy(IE)
        endif
        PT(1,IG)=P(1,IE)
        PT(2,IG)=P(2,IE)
        PT(3,IG)=P(3,IE)
        ET(IG)=E(IE)
        LT(IG)=LB(IE)
        PROT(IG)=PROPER(IE)
        dptemp(IG)=dpertp(IE)
        ELSE
        I0=IC-MASS
        RT(1,IG)=RPION(1,I0,IRUN)
        RT(2,IG)=RPION(2,I0,IRUN)
        RT(3,IG)=RPION(3,I0,IRUN)
        if(nt.eq.ntmax) then
           fttemp(IG)=ftpisv(I0,IRUN)
           tft(IG)=tfdpi(I0,IRUN)
        endif
        PT(1,IG)=PPION(1,I0,IRUN)
        PT(2,IG)=PPION(2,I0,IRUN)
        PT(3,IG)=PPION(3,I0,IRUN)
        ET(IG)=EPION(I0,IRUN)
        LT(IG)=LPION(I0,IRUN)
        PROT(IG)=PROPI(I0,IRUN)
        dptemp(IG)=dppion(I0,IRUN)
        ENDIF
10001   CONTINUE
        IL=0
        DO 10003 IRUN=1,NUM
        MASSR(IRUN)=MASSRN(IRUN)
        IL=IL+MASSR(IRUN-1)
        DO 10002 IM=1,MASSR(IRUN)
        IN=IL+IM
        R(1,IN)=RT(1,IN)
        R(2,IN)=RT(2,IN)
        R(3,IN)=RT(3,IN)
        if(nt.eq.ntmax) then
           ftsv(IN)=fttemp(IN)
           tfdcy(IN)=tft(IN)
        endif
        P(1,IN)=PT(1,IN)
        P(2,IN)=PT(2,IN)
        P(3,IN)=PT(3,IN)
        E(IN)=ET(IN)
        LB(IN)=LT(IN)
        PROPER(IN)=PROT(IN)
        dpertp(IN)=dptemp(IN)
       IF(LB(IN).LT.1.OR.LB(IN).GT.2)ID(IN)=0
10002   CONTINUE
10003 CONTINUE
      RETURN
      END
