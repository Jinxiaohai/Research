      subroutine addhad
      PARAMETER (MAXSTR=150001,MAXR=1,xmd=1.8756)
      double precision  smearp,smearh
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)
      COMMON /ARPRC/ ITYPAR(MAXSTR),
     &     GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &     PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &     XMAR(MAXSTR)
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
      COMMON /smearz/smearp,smearh
      COMMON/RNDF77/NSEED
      common /para8/ idpert,npertd,idxsec
      SAVE   
      np0=IAINT2(1)
      DO i=1,np0
         dpertp(I)=1.
      ENDDO
      nadd=0
      tau0=ARPAR1(1)
      DO 100 i=np0+1,np0+nadd
         ITYPAR(I)=42
         dpertp(I)=1./float(nadd)
         GXAR(I)=5.*(1.-2.*RANART(NSEED))
         GYAR(I)=5.*(1.-2.*RANART(NSEED))
         GZAR(I)=2.*(1.-2.*RANART(NSEED))
         FTAR(I)=0.
         PXAR(I)=1.
         PYAR(I)=0.
         PZAR(I)=1.
         XMAR(I)=xmd
         PEAR(I)=sqrt(PXAR(I)**2+PYAR(I)**2+PZAR(I)**2+XMAR(I)**2)
         RAP=asinh(PZAR(I)/sqrt(XMAR(I)**2+PXAR(I)**2+PYAR(I)**2))
         VX=PXAR(I)/PEAR(I)
         VY=PYAR(I)/PEAR(I)
         TAUI=FTAR(I)+TAU0
         FTAR(I)=TAUI*COSH(RAP)
         GXAR(I)=GXAR(I)+VX*TAU0*COSH(RAP)
         GYAR(I)=GYAR(I)+VY*TAU0*COSH(RAP)
         GZAR(I)=TAUI*SINH(RAP)+GZAR(I)
         zsmear=sngl(smearh)*(2.*RANART(NSEED)-1.)
         GZAR(I)=GZAR(I)+zsmear
 100  CONTINUE
      IAINT2(1)=IAINT2(1)+nadd
      if(nadd.ge.1.and.idpert.ne.1.and.idpert.ne.2) then
         write(16,*) 'IDPERT must be 1 or 2 to add initial hadrons,
     1 set NPERTD to 0 if you do not need perturbative deuterons'
         stop
      endif
      if(IAINT2(1).gt.MAXSTR) then
         write(16,*) 'Too many initial hadrons, array size is exceeded!'
         stop
      endif
      return
      end
