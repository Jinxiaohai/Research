      SUBROUTINE HJAN2B
      PARAMETER (MAXPTN=400001)
      PARAMETER (MAXSTR=150001)
      PARAMETER (DR = 0.2, DT = 0.2)
      DIMENSION DNRG2B(50), dtg2b(-24:25)
      DIMENSION SNRG2B(50), stg2b(-24:25)
      DOUBLE PRECISION  GX5, GY5, GZ5, FT5, PX5, PY5, PZ5, E5, XMASS5
      DOUBLE PRECISION  ATAUI, ZT1, ZT2, ZT3
      COMMON /PARA1/ MUL
      COMMON /prec2/GX5(MAXPTN),GY5(MAXPTN),GZ5(MAXPTN),FT5(MAXPTN),
     &   PX5(MAXPTN), PY5(MAXPTN), PZ5(MAXPTN), E5(MAXPTN),
     &   XMASS5(MAXPTN), ITYP5(MAXPTN)
      COMMON /ilist8/ LSTRG1(MAXPTN), LPART1(MAXPTN)
      COMMON /SREC1/ NSP, NST, NSI
      COMMON /SREC2/ATAUI(MAXSTR),ZT1(MAXSTR),ZT2(MAXSTR),ZT3(MAXSTR)
      COMMON/hjcrdn/YP(3,300),YT(3,300)
      COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &   K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &   PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
      COMMON /AREVT/ IAEVT, IARUN, MISS
      COMMON /AROUT/ IOUT
      SAVE   
      DATA IW/0/
      IF (isevt .EQ. IAEVT .AND. isrun .EQ. IARUN) THEN
         DO 1001 I = 1, 50
            DNRG2B(I) = SNRG2B(I)
            dtg2b(I - 25) = stg2b(I - 25)
 1001    CONTINUE
      ELSE
         DO 1002 I = 1, 50
            SNRG2B(I) = DNRG2B(I)
            stg2b(I - 25) = dtg2b(I - 25)
 1002    CONTINUE
         isevt = IAEVT
         isrun = IARUN
         IW = IW + 1
      END IF
      DO 1003 I = 1, MUL
         J = LSTRG1(I)
         GX0 = sngl(ZT1(J))
         GY0 = sngl(ZT2(J))
         R0 = SQRT((sngl(GX5(I)) - GX0)**2 + (sngl(GY5(I)) - GY0)**2)
         IR = 1 + int(R0 / DR)
         IF (IR .GT. 50 .or. IR .LT. 1) GOTO 100
         DNRG2B(IR) = DNRG2B(IR) + 1.0
 100     CONTINUE
         diff2=sngl(FT5(I)**2 - GZ5(I)**2)
         if(diff2.lt.0.) then
            write(6,*) '4:I,ft5,gz5,diff2=',I,ft5(i),gz5(i),diff2
            TAU7=1e-6
         else
            TAU7 = SQRT(diff2)
         endif
         DTAU=TAU7 - sngl(ATAUI(J))
         IT = 1 + int(DTAU / DT)
         IF (IT .GT. 25 .OR. IT .LT. -24) GOTO 200
         dtg2b(IT) = dtg2b(IT) + 1.0
 200     CONTINUE
 1003 CONTINUE
      RETURN
      END
