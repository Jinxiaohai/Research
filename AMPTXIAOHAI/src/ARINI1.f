      SUBROUTINE ARINI1
      PARAMETER (MAXSTR=150001)
      double precision  smearp,smearh
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)
      COMMON /ARPRC/ ITYPAR(MAXSTR),
     &     GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &     PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &     XMAR(MAXSTR)
      COMMON /smearz/smearp,smearh
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      common/anim/nevent,isoft,isflag,izpc
      common /nzpc/nattzp
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      COMMON/RNDF77/NSEED
      common /para8/ idpert,npertd,idxsec
      SAVE   
      OPEN (91, FILE = 'ana/deuteron_processes.dat', 
     1     STATUS = 'UNKNOWN')
      if(idpert.eq.1.or.idpert.eq.2) then
         OPEN (90, FILE = 'ana/ampt_pert.dat', STATUS = 'UNKNOWN')
      endif
      TAU0 = ARPAR1(1)
      NP = IAINT2(1)
      if(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) then
         if(NP.gt.nattzp) then
         do 1001 I = nattzp+1, NP
            if((XMAR(I)**2+PXAR(I)**2+PYAR(I)**2).gt.0.) then
               RAP=asinh(PZAR(I)/sqrt(XMAR(I)**2+PXAR(I)**2+PYAR(I)**2))
            else
               PRINT *, ' IN ARINI1 mt=0'
               RAP = 1000000.0*sign(1.,PZAR(I))
            endif
            VX = PXAR(I) / PEAR(I)
            VY = PYAR(I) / PEAR(I)
            FTAR(I) = TAU0 * COSH(RAP)
            GXAR(I) = GXAR(I) + VX * FTAR(I)
            GYAR(I) = GYAR(I) + VY * FTAR(I)
            GZAR(I) = TAU0 * SINH(RAP)
            if(PXAR(I).eq.0.and.PYAR(I).eq.0
     2           .and.(ITYPAR(I).eq.2112.or.ITYPAR(I).eq.2212)) then
              if((PEAR(I)/HINT1(6).gt.0.99.and.PEAR(I)/HINT1(6).lt.1.01)
     1 .or.(PEAR(I)/HINT1(7).gt.0.99.and.PEAR(I)/HINT1(7).lt.1.01)) then
               TAUI=1.E-20
               FTAR(I)=TAUI*COSH(RAP)
               GZAR(I)=TAUI*SINH(RAP)
               endif
            endif
 1001    continue
         endif
      else
         DO 1002 I = 1, NP
            if((XMAR(I)**2+PXAR(I)**2+PYAR(I)**2).gt.0.) then
               RAP=asinh(PZAR(I)/sqrt(XMAR(I)**2+PXAR(I)**2+PYAR(I)**2))
            else
               PRINT *, ' IN ARINI1 mt=0'
               RAP = 1000000.0*sign(1.,PZAR(I))
            endif
            VX = PXAR(I) / PEAR(I)
            VY = PYAR(I) / PEAR(I)
            TAUI = FTAR(I) + TAU0
            FTAR(I) = TAUI * COSH(RAP)
            GXAR(I) = GXAR(I) + VX * TAU0 * COSH(RAP)
            GYAR(I) = GYAR(I) + VY * TAU0 * COSH(RAP)
            GZAR(I) = TAUI * SINH(RAP)
            zsmear=sngl(smearh)*(2.*RANART(NSEED)-1.)
            GZAR(I)=GZAR(I)+zsmear
            if(PXAR(I).eq.0.and.PYAR(I).eq.0
     2           .and.(ITYPAR(I).eq.2112.or.ITYPAR(I).eq.2212)) then
              if((PEAR(I)/HINT1(6).gt.0.99.and.PEAR(I)/HINT1(6).lt.1.01)
     1 .or.(PEAR(I)/HINT1(7).gt.0.99.and.PEAR(I)/HINT1(7).lt.1.01)) then
               TAUI=1.E-20
               FTAR(I)=TAUI*COSH(RAP)
               GZAR(I)=TAUI*SINH(RAP)+zsmear
               endif
            endif
 1002    CONTINUE
      endif
      call addhad
      RETURN
      END
