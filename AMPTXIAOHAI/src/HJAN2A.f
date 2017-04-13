      SUBROUTINE HJAN2A
      PARAMETER (DGX = 0.2, DGY = 0.2, DT = 0.2)
      PARAMETER (MAXPTN=400001,MAXSTR=150001)
      DIMENSION dgxp2a(50), dgyp2a(50), dtp2a(50)
      DIMENSION dgxg2a(50), dgyg2a(50), dtg2a(50)
      DIMENSION sgxp2a(50), sgyp2a(50), stp2a(50)
      DIMENSION sgxg2a(50), sgyg2a(50), stg2a(50)
      DOUBLE PRECISION  GX5, GY5, GZ5, FT5, PX5, PY5, PZ5, E5, XMASS5
      COMMON /PARA1/ MUL
      COMMON /prec2/GX5(MAXPTN),GY5(MAXPTN),GZ5(MAXPTN),FT5(MAXPTN),
     &   PX5(MAXPTN), PY5(MAXPTN), PZ5(MAXPTN), E5(MAXPTN),
     &   XMASS5(MAXPTN), ITYP5(MAXPTN)
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      COMMON/hjcrdn/YP(3,300),YT(3,300)
      COMMON/HJJET1/NPJ(300),KFPJ(300,500),PJPX(300,500),
     &   PJPY(300,500),PJPZ(300,500),PJPE(300,500),
     &   PJPM(300,500),NTJ(300),KFTJ(300,500),
     &   PJTX(300,500),PJTY(300,500),PJTZ(300,500),
     &   PJTE(300,500),PJTM(300,500)
      COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &   K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &   PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
      COMMON /AREVT/ IAEVT, IARUN, MISS
      COMMON /AROUT/ IOUT
      SAVE   
      DATA IW/0/
      IF (isevt .EQ. IAEVT .AND. isrun .EQ. IARUN) THEN
         DO 1001 I = 1, 50
            dgxp2a(I) = sgxp2a(I)
            dgyp2a(I) = sgyp2a(I)
            dtp2a(I) = stp2a(I)
            dgxg2a(I) = sgxg2a(I)
            dgyg2a(I) = sgyg2a(I)
            dtg2a(I) = stg2a(I)
 1001    CONTINUE
      ELSE
         DO 1002 I = 1, 50
            sgxp2a(I) = dgxp2a(I)
            sgyp2a(I) = dgyp2a(I)
            stp2a(I) = dtp2a(I)
            sgxg2a(I) = dgxg2a(I)
            sgyg2a(I) = dgyg2a(I)
            stg2a(I) = dtg2a(I)
 1002    CONTINUE
         isevt = IAEVT
         isrun = IARUN
         IW = IW + 1
      END IF
      DO 1004 I = 1, IHNT2(1)
         DO 1003 J = 1, NPJ(I)
            IF (KFPJ(I, J) .NE. 21) THEN
               IGX = 1 + int(ABS(YP(1, I)) / DGX)
               IF (IGX .GT. 50 .or. IGX .LT. 1) GOTO 100
               dgxp2a(IGX) = dgxp2a(IGX) + 1.0
 100           CONTINUE
               IGY = 1 + int(ABS(YP(2, I)) / DGY)
               IF (IGY .GT. 50 .or. IGY .LT. 1) GOTO 200
               dgyp2a(IGY) = dgyp2a(IGY) + 1.0
 200           CONTINUE
               IT = 1
               dtp2a(IT) = dtp2a(IT) + 1.0
            END IF
 1003    CONTINUE
 1004 CONTINUE
      DO 1006 I = 1, IHNT2(3)
         DO 1005 J = 1, NTJ(I)
            IF (KFTJ(I, J) .NE. 21) THEN
               IGX = 1 + int(ABS(YT(1, I)) / DGX)
               IF (IGX .GT. 50 .or. IGX .LT. 1) GOTO 300
               dgxp2a(IGX) = dgxp2a(IGX) + 1.0
 300           CONTINUE
               IGY = 1 + int(ABS(YT(2, I)) / DGY)
               IF (IGY .GT. 50 .or. IGY .LT. 1) GOTO 400
               dgyp2a(IGY) = dgyp2a(IGY) + 1.0
 400           CONTINUE
               IT = 1
               dtp2a(IT) = dtp2a(IT) + 1.0
            END IF
 1005    CONTINUE
 1006 CONTINUE
      DO 1008 I = 1, NSG
         DO 1007 J = 1, NJSG(I)
            IF (K2SG(I, J) .NE. 21) THEN
               IGX = 1 + int(ABS(0.5 * 
     &            (YP(1, IASG(I, 1)) + YT(1, IASG(I, 2)))) / DGX)
               IF (IGX .GT. 50 .or. IGX .LT. 1) GOTO 500
               dgxp2a(IGX) = dgxp2a(IGX) + 1.0
 500           CONTINUE
               IGY = 1 + int(ABS(0.5 * 
     &            (YP(2, IASG(I, 1)) + YT(2, IASG(I, 2)))) / DGY)
               IF (IGY .GT. 50 .or. IGY .LT. 1) GOTO 600
               dgyp2a(IGY) = dgyp2a(IGY) + 1.0
 600           CONTINUE
               IT = 1
               dtp2a(IT) = dtp2a(IT) + 1.0               
            END IF
 1007    CONTINUE
 1008 CONTINUE
      DO 1009 I = 1, MUL
         IGX = 1 + int(ABS(sngl(GX5(I))) / DGX)
         IF (IGX .GT. 50 .or. IGX .LT. 1) GOTO 700
         dgxg2a(IGX) = dgxg2a(IGX) + 1.0
         dgxp2a(IGX) = dgxp2a(IGX) + 1.0
 700     CONTINUE
         IGY = 1 + int(ABS(sngl(GY5(I))) / DGY)
         IF (IGY .GT. 50 .or. IGY .LT. 1) GOTO 800
         dgyg2a(IGY) = dgyg2a(IGY) + 1.0
         dgyp2a(IGY) = dgyp2a(IGY) + 1.0
 800     CONTINUE
         diff2=sngl(FT5(I)**2 - GZ5(I)**2)
         if(diff2.lt.0.) then
            write(6,*) '3:I,ft5,gz5,diff2=',I,ft5(i),gz5(i),diff2
            IT=1
         else
            IT = 1 + int(SQRT(diff2)/DT)
         endif
         IF (IT .GT. 50 .or. IT .LT. 1) GOTO 900
         dtg2a(IT) = dtg2a(IT) + 1.0
         dtp2a(IT) = dtp2a(IT) + 1.0
 900     CONTINUE
 1009 CONTINUE
      RETURN
      END
