      PROGRAM AMPT
      double precision xmp, xmu, alpha, rscut2, cutof2, dshadow
      double precision smearp,smearh,dpcoal,drcoal,ecritl
      CHARACTER FRAME*8, PROJ*8, TARG*8
      character*25 amptvn
      COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
      COMMON /HPARNT/HIPR1(100), IHPR2(50), HINT1(100), IHNT2(50)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)
      COMMON /AROUT/ IOUT
      COMMON /AREVT/ IAEVT, IARUN, MISS
      COMMON /smearz/smearp,smearh
      COMMON/RNDF77/NSEED
      common/anim/nevent,isoft,isflag,izpc
      common /coal/dpcoal,drcoal,ecritl
      common/snn/efrm,npart1,npart2,epsiPz,epsiPt,PZPROJ,PZTARG
      common /para2/ xmp, xmu, alpha, rscut2, cutof2
      common /para7/ ioscar,nsmbbbar,nsmmeson
      common /para8/ idpert,npertd,idxsec
      common /rndm3/ iseedp
      COMMON /RUN/ NUM
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &   IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
      common/oscar1/iap,izp,iat,izt
      common/oscar2/FRAME,amptvn
      common/resdcy/NSAV,iksdcy
      common/phidcy/iphidcy,pttrig,ntrig,maxmiss,ipi0dcy
      common/embed/iembed,nsembd,pxqembd,pyqembd,xembd,yembd,
     1     psembd,tmaxembd,phidecomp
      common/cmsflag/dshadow,ishadow
      common /phiHJ/iphirp,phiRP
      EXTERNAL HIDATA, PYDATA, LUDATA, ARDATA, PPBDAT, zpcbdt
      SAVE   
      OPEN (24, FILE = 'input.ampt', STATUS = 'UNKNOWN')
      OPEN (12, FILE = 'ana/version', STATUS = 'UNKNOWN')
      READ (24, *) EFRM
      READ (24, 111) FRAME
      READ (24, 111) PROJ
      READ (24, 111) TARG
      READ (24, *) IAP
      READ (24, *) IZP
      READ (24, *) IAT
      READ (24, *) IZT
      READ (24, *) NEVNT
      READ (24, *) BMIN
      READ (24, *) BMAX
      READ (24, *) isoft
      READ (24, *) NTMAX
      READ (24, *) DT
      READ (24, *) PARJ(41)
      READ (24, *) PARJ(42)
      READ (24, *) ipop
      if(ipop.eq.1) IHPR2(11)=3
      READ (24, *) PARJ(5)
      READ (24, *) IHPR2(6)
      READ (24, *) IHPR2(4)
      READ (24, *) HIPR1(14)
      READ (24, *) HIPR1(8)
      READ (24, *) xmu
      READ (24, *) izpc
      READ (24, *) alpha
      READ (24, *) dpcoal
      READ (24, *) drcoal
      READ (24, *) ihjsed
      READ (24, *) nseed
      READ (24, *) iseedp
      READ (24, *) iksdcy
      READ (24, *) iphidcy
      READ (24, *) ipi0dcy
      READ (24, *) ioscar
      READ (24, *) idpert
      READ (24, *) npertd
      READ (24, *) idxsec
      READ (24, *) pttrig
      READ (24, *) maxmiss
      READ (24, *) IHPR2(2)
      READ (24, *) IHPR2(5)
      READ (24, *) iembed
      READ (24, *) pxqembd, pyqembd
      READ (24, *) xembd, yembd
      READ (24, *) nsembd,psembd,tmaxembd
      READ (24, *) ishadow
      READ (24, *) dshadow
      READ (24, *) iphirp
      CLOSE (24)
 111  format(a8)
      if(isoft.eq.1) then
         amptvn = '1.26t7 (Default)'
      elseif(isoft.eq.4) then
         amptvn = '2.26t7 (StringMelting)'
      else
         amptvn = 'Test-Only'
      endif
      WRITE(6,50) amptvn
      WRITE(12,50) amptvn
 50   FORMAT(' '/
     &11X,'##################################################'/1X,
     &10X,'#      AMPT (A Multi-Phase Transport) model      #'/1X,
     &10X,'#          Version ',a25,                  '     #'/1X,
     &10X,'#               10/28/2016                       #'/1X,
     &10X,'##################################################'/1X,
     &10X,' ')
      if(ihjsed.eq.11) then
         PRINT *,
     1 '# Read in NSEED in HIJING at run time (e.g. 20030819):'
      endif
      READ (*, *) nseedr
      if(ihjsed.eq.11) then
         nseed=nseedr
      endif
      if(ihjsed.eq.11) then      
         PRINT *, '#   read in: ', nseed
         WRITE(12,*) '# Read in NSEED in HIJING at run time:',nseed
      endif
      CLOSE(12)
      NSEED=2*NSEED+1
      CALL SRAND(NSEED)
      IHPR2(10) = 1
      ARPAR1(1) = 0.7
      smearp=0d0
      IAmax=max(iap,iat)
      smearh=1.2d0*IAmax**0.3333d0/(dble(EFRM)/2/0.938d0)
      nevent=NEVNT
      OPEN (16, FILE = 'ana/ampt.dat', STATUS = 'UNKNOWN')
      OPEN (14, FILE = 'ana/zpc.dat', STATUS = 'UNKNOWN')
      CALL HIJSET(EFRM, FRAME, PROJ, TARG, IAP, IZP, IAT, IZT)
      CALL ARTSET
      CALL INIZPC
       DO 2000 J = 1, NEVNT
          IAEVT = J
          DO 1000 K = 1, NUM
             IARUN = K
             IF (IAEVT .EQ. NEVNT .AND. IARUN .EQ. NUM) THEN
                IOUT = 1
             END IF
             PRINT *, ' EVENT ', J, ', RUN ', K
             imiss=0
 100         CALL HIJING(FRAME, BMIN, BMAX)
             IAINT2(1) = NATT             
           if(J.eq.-2) then 
              write(98,*) HIPR1
              write(98,*) ' '
              write(98,*) IHPR2
              write(98,*) ' '
              write(98,*) (HINT1(i),i=1,20)
              write(98,*) ' '
              write(98,*) (HINT1(i),i=21,40)
              write(98,*) ' '
              write(98,*) (HINT1(i),i=41,60)
              write(98,*) ' '
              write(98,*) (HINT1(i),i=61,80)
              write(98,*) ' '
              write(98,*) (HINT1(i),i=81,100)
              write(98,*) ' '
              write(98,*) IHNT2
           endif
             call getnp
             IF (IHPR2(20) .EQ. 0) GOTO 2000
             if(NATT.eq.0) then
                imiss=imiss+1
                if(imiss.le.20) then
                   write(6,*) 'repeated event: natt=0,j,imiss=',j,imiss
                   goto 100
                else
                   write(6,*) 'missed event: natt=0,j=',j
                   goto 2000
                endif
             endif
             CALL ARINI
             CALL ARINI2(K)
 1000     CONTINUE
          CALL ARTAN1
          CALL ARTMN
          CALL ARTAN2
 2000  CONTINUE
       CALL ARTOUT(NEVNT)
       STOP
       END
