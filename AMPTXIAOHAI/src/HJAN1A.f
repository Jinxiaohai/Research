      SUBROUTINE HJAN1A
      PARAMETER (MAXPTN=400001)
      PARAMETER (DGX = 0.2, DGY = 0.2, DT = 0.2)
      DIMENSION dgxg1a(50), dgyg1a(50), dtg1a(50)
      DIMENSION sgxg1a(50), sgyg1a(50), stg1a(50)
      DOUBLE PRECISION  GX5, GY5, GZ5, FT5, PX5, PY5, PZ5, E5, XMASS5
      COMMON /PARA1/ MUL
      COMMON /prec2/GX5(MAXPTN),GY5(MAXPTN),GZ5(MAXPTN),FT5(MAXPTN),
     &   PX5(MAXPTN), PY5(MAXPTN), PZ5(MAXPTN), E5(MAXPTN),
     &   XMASS5(MAXPTN), ITYP5(MAXPTN)
      COMMON /AREVT/ IAEVT, IARUN, MISS
      COMMON /AROUT/ IOUT
      SAVE   
      DATA IW/0/
      IF (isevt .EQ. IAEVT .AND. isrun .EQ. IARUN) THEN
         DO 1001 I = 1, 50
            dgxg1a(I) = sgxg1a(I)
            dgyg1a(I) = sgyg1a(I)
            dtg1a(I) = stg1a(I)
 1001    CONTINUE
      ELSE
         DO 1002 I = 1, 50
            sgxg1a(I) = dgxg1a(I)
            sgyg1a(I) = dgyg1a(I)
            stg1a(I) = dtg1a(I)
 1002    CONTINUE
         isevt = IAEVT
         isrun = IARUN
         IW = IW + 1
      END IF
      DO 1003 I = 1, MUL
         IGX = 1 + int(sngl(ABS(GX5(I))) / DGX)
         IF (IGX .GT. 50 .or. IGX .LT. 1) GOTO 100
         dgxg1a(IGX) = dgxg1a(IGX) + 1.0
 100     CONTINUE
         IGY = 1 + int(sngl(ABS(GY5(I))) / DGY)
         IF (IGY .GT. 50 .or. IGY .LT. 1) GOTO 200
         dgyg1a(IGY) = dgyg1a(IGY) + 1.0
 200     CONTINUE
         diff2=sngl(FT5(I)**2 - GZ5(I)**2)
         if(diff2.lt.0.) then
            write(6,*) '1:I,ft5,gz5,diff2=',I,ft5(i),gz5(i),diff2
            IT=1
         else
            IT = 1 + int(SQRT(diff2)/DT)
         endif
         IF (IT .GT. 50 .or. IT .LT. 1) GOTO 300
         dtg1a(IT) = dtg1a(IT) + 1.0
 300     CONTINUE
 1003 CONTINUE
      CALL HJAN1B
      RETURN
      END
