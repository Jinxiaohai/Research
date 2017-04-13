      SUBROUTINE ARTSET
      PARAMETER (AMU= 0.9383,nxymax=10001)
      double precision dpcoal,drcoal,ecritl
      INTEGER ZTA, ZPR
      common  /gg/      dx,dy,dz,dpx,dpy,dpz
      common  /zz/      zta,zpr
      COMMON  /RUN/     NUM
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &   IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
      COMMON /INPUT3/ PLAB, ELAB, ZEROPT, B0, BI, BM, DENCUT, CYCBOX
      common /imulst/ iperts
      common /coal/dpcoal,drcoal,ecritl
      common/anim/nevent,isoft,isflag,izpc
      common /para7/ ioscar,nsmbbbar,nsmmeson
      common/embed/iembed,nsembd,pxqembd,pyqembd,xembd,yembd,
     1     psembd,tmaxembd,phidecomp
      common/xyembed/nxyjet,xyjet(nxymax,2)
      SAVE   
      ecritl=1.d0
      MASSTA=1
      ZTA=1
      MASSPR=1
      ZPR=1
      PLAB=14.6 
      IPLAB=2
      if(iplab.eq.2)then
         elab=sqrt(plab**2+amu**2)-amu
      else
         elab=plab
      endif
      elab=elab*1000.
      ZEROPT=0.
      ISEED=700721
      iperts=0
      MANYB=1
      B0=1
      BI=0
      BM=0
      ICOLL=-1
      NUM=1
      INSYS=1
      IPOT=3
      MODE=0
      IF(ICOLL.EQ.-1)IPOT=0
      DX=2.73
      DY=2.73
      DZ=2.73
      DPX=0.6
      DPY=0.6
      DPZ=0.6
      IAVOID=1
      IMOMEN=1
      if(icoll.eq.-1)imomen=3
      NFREQ=10
      ICFLOW=0
      ICRHO=0
      ICOU=0
      KPOTEN=0
      KMUL=1
      DENCUT=15
      CYCBOX=0
      if(ioscar.eq.2.or.ioscar.eq.3) then
         OPEN (92,FILE='ana/parton-initial-afterPropagation.dat',
     1        STATUS = 'UNKNOWN')
      endif
      if(ioscar.eq.3) then
         OPEN (95,FILE='ana/parton-collisionsHistory.dat',
     1        STATUS='UNKNOWN')
         OPEN (96,FILE='ana/minijet-initial-beforePropagation.dat',
     1        STATUS='UNKNOWN')
         if(isoft.eq.4.or.isoft.eq.5) then
            OPEN (85,FILE='ana/parton-after-coalescence.dat',
     1           STATUS='UNKNOWN')
         endif
      endif
      OPEN (94,FILE='ana/npart-xy.dat',STATUS='UNKNOWN')
      if(iembed.eq.3.or.iembed.eq.4) then
         OPEN (97,FILE='embed-jet-xy.txt',STATUS = 'UNKNOWN')
         read(97,*) nxyjet
         if(nevent.gt.nxyjet) then
            if(nxyjet.gt.nxymax) then
               print *, 'Too many lines in embed-jet-xy.txt: 
     1 increase value of the parameter nxymax'
               stop
            elseif(nxyjet.le.0) then
               print *, 'Check number of entries in embed-jet-xy.txt' 
               stop
            endif
            do ixy=1,nxyjet
               read(97,*) xyjet(ixy,1),xyjet(ixy,2)
            enddo
         endif
      endif
      RETURN
      END
