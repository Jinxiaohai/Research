      SUBROUTINE ARTORD
      PARAMETER (MAXSTR=150001,MAXR=1)
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)
      COMMON /ARPRC/ ITYPAR(MAXSTR),
     &     GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &     PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &     XMAR(MAXSTR)
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
      DIMENSION dptemp(MAXSTR)
      DIMENSION ITYP0(MAXSTR), 
     &   GX0(MAXSTR), GY0(MAXSTR), GZ0(MAXSTR), FT0(MAXSTR),
     &   PX0(MAXSTR), PY0(MAXSTR), PZ0(MAXSTR), EE0(MAXSTR),
     &   XM0(MAXSTR)
      DIMENSION INDX(MAXSTR)
      EXTERNAL ARINDX
      SAVE   
      NPAR = 0
      NP = IAINT2(1)
      DO 1001 I = 1, NP
         ITYP0(I) = ITYPAR(I)
         GX0(I) = GXAR(I)
         GY0(I) = GYAR(I)
         GZ0(I) = GZAR(I)
         FT0(I) = FTAR(I)
         PX0(I) = PXAR(I)
         PY0(I) = PYAR(I)
         PZ0(I) = PZAR(I)
         EE0(I) = PEAR(I)
         XM0(I) = XMAR(I)
         dptemp(I) = dpertp(I)
 1001 CONTINUE
      CALL ARINDX(MAXSTR, NP, FT0, INDX)
      DO 1002 I = 1, NP
         NPAR = NPAR + 1
         ITYPAR(NPAR) = ITYP0(INDX(I))
         GXAR(NPAR) = GX0(INDX(I))
         GYAR(NPAR) = GY0(INDX(I))
         GZAR(NPAR) = GZ0(INDX(I))
         FTAR(NPAR) = FT0(INDX(I))
         PXAR(NPAR) = PX0(INDX(I))
         PYAR(NPAR) = PY0(INDX(I))
         PZAR(NPAR) = PZ0(INDX(I))
         PEAR(NPAR) = EE0(INDX(I))
         XMAR(NPAR) = XM0(INDX(I))
         dpertp(NPAR)=dptemp(INDX(I))
 1002 CONTINUE
      IAINT2(1) = NPAR
      RETURN
      END
