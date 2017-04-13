      SUBROUTINE HJAN1B
      PARAMETER (MAXPTN=400001,MAXSTR=150001)
      PARAMETER (DR = 0.2, DT = 0.2)
      DIMENSION DNRG1B(50), dtg1b(50)
      DIMENSION SNRG1B(50), stg1b(50)
      DOUBLE PRECISION  GX5, GY5, GZ5, FT5, PX5, PY5, PZ5, E5, XMASS5
      COMMON /PARA1/ MUL
      COMMON /prec2/GX5(MAXPTN),GY5(MAXPTN),GZ5(MAXPTN),FT5(MAXPTN),
     &   PX5(MAXPTN), PY5(MAXPTN), PZ5(MAXPTN), E5(MAXPTN),
     &   XMASS5(MAXPTN), ITYP5(MAXPTN)
      COMMON /ilist8/ LSTRG1(MAXPTN), LPART1(MAXPTN)
      COMMON /SREC1/ NSP, NST, NSI
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
            DNRG1B(I) = SNRG1B(I)
            dtg1b(I) = stg1b(I)
 1001    CONTINUE
      ELSE
         DO 1002 I = 1, 50
            SNRG1B(I) = DNRG1B(I)
            stg1b(I) = dtg1b(I)
 1002    CONTINUE
         isevt = IAEVT
         isrun = IARUN
         IW = IW + 1
      END IF
      DO 1003 I = 1, MUL
         J = LSTRG1(I)
         IF (J .LE. NSP) THEN
            K = J
            GX0 = YP(1, J)
            GY0 = YP(2, J)
         ELSE IF (J .LE. NSP + NST) THEN
            K = J - NSP
            GX0 = YT(1, K)
            GY0 = YT(2, K)
         ELSE
            K = J - NSP - NST
            GX0 = 0.5 * (YP(1, IASG(K, 1)) + YT(1, IASG(K, 2)))
            GY0 = 0.5 * (YP(2, IASG(K, 1)) + YT(2, IASG(K, 2)))
         END IF
         R0 = SQRT((sngl(GX5(I)) - GX0)**2 + (sngl(GY5(I)) - GY0)**2)
         IR = 1 + int(R0 / DR)
         IF (IR .GT. 50 .or. IR .LT. 1) GOTO 100
         DNRG1B(IR) = DNRG1B(IR) + 1.0
 100     CONTINUE
         diff2=sngl(FT5(I)**2 - GZ5(I)**2)
         if(diff2.lt.0.) then
            write(6,*) '5:I,ft5,gz5,diff2=',I,ft5(i),gz5(i),diff2
            TAU7 = 1e-6
         else
            TAU7 = SQRT(diff2)
         endif
         IT = 1 + int(TAU7 / DT)
         IF (IT .GT. 50 .or. IT .LT. 1) GOTO 200
         dtg1b(IT) = dtg1b(IT) + 1.0
 200     CONTINUE
 1003 CONTINUE
      RETURN
      END
