      SUBROUTINE HJANA1
      PARAMETER (YMAX = 1.0, YMIN = -1.0)
      PARAMETER (DMT = 0.05, DY = 0.2)
      PARAMETER (DR = 0.2)
      PARAMETER (MAXPTN=400001,MAXSTR=150001)
      DIMENSION dyp1(50), DMYP1(200), DEYP1(50)
      DIMENSION dyg1(50), DMYG1(200), DEYG1(50)
      DIMENSION SNYP1(50), SMYP1(200), SEYP1(50)
      DIMENSION SNYG1(50), SMYG1(200), SEYG1(50)
      DIMENSION dnrpj1(50), dnrtg1(50), dnrin1(50),
     &   dnrtt1(50)
      DIMENSION dyg1c(50), dmyg1c(50), deyg1c(50)
      DIMENSION snrpj1(50), snrtg1(50), snrin1(50),
     &   snrtt1(50)
      DIMENSION snyg1c(50), smyg1c(50), seyg1c(50)
      DOUBLE PRECISION  GX0, GY0, GZ0, FT0, PX0, PY0, PZ0, E0, XMASS0
      COMMON /PARA1/ MUL
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
      COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &     PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &     XMASS0(MAXPTN), ITYP0(MAXPTN)
      COMMON /AREVT/ IAEVT, IARUN, MISS
      COMMON /AROUT/ IOUT
      SAVE   
      DATA IW/0/
      IF (isevt .EQ. IAEVT .AND. isrun .EQ. IARUN) THEN
         DO 1001 I = 1, 200
            DMYP1(I) = SMYP1(I)
            DMYG1(I) = SMYG1(I)
 1001    CONTINUE
         DO 1002 I = 1, 50
            dyp1(I) = SNYP1(I)
            DEYP1(I) = SEYP1(I)
            dyg1(I) = SNYG1(I)
            DEYG1(I) = SEYG1(I)
            dnrpj1(I) = snrpj1(I)
            dnrtg1(I) = snrtg1(I)
            dnrin1(I) = snrin1(I)
            dnrtt1(I) = snrtt1(I)
            dyg1c(I) = snyg1c(I)
            dmyg1c(I) = smyg1c(I)
            deyg1c(I) = seyg1c(I)
 1002    CONTINUE
         nsubp = nsubpS
         nsubg = nsubgS
         NISG = NISGS
      ELSE
         DO 1003 I = 1, 200
            SMYP1(I) = DMYP1(I)
            SMYG1(I) = DMYG1(I)
 1003    CONTINUE
         DO 1004 I = 1, 50
            SNYP1(I) = dyp1(I)
            SEYP1(I) = DEYP1(I)
            SNYG1(I) = dyg1(I)
            SEYG1(I) = DEYG1(I)
            snrpj1(I) = dnrpj1(I)
            snrtg1(I) = dnrtg1(I)
            snrin1(I) = dnrin1(I)
            snrtt1(I) = dnrtt1(I)
            snyg1c(I) = dyg1c(I)
            smyg1c(I) = dmyg1c(I)
            seyg1c(I) = deyg1c(I)
 1004    CONTINUE
         nsubpS = nsubp
         nsubgS = nsubg
         NISGS = NISG
         isevt = IAEVT
         isrun = IARUN
         IW = IW + 1
      END IF
      DO 1006 I = 1, IHNT2(1)
         DO 1005 J = 1, NPJ(I)
            ITYP = KFPJ(I, J)
            PX = PJPX(I, J)
            PY = PJPY(I, J)
            PZ = PJPZ(I, J)
            PE = PJPE(I, J)
            PM = PJPM(I, J)
            XMT = SQRT(PX ** 2 + PY ** 2 + PM ** 2)
            DXMT = XMT - PM
            if(XMT.gt.0.) then
               RAP=asinh(PZ/XMT)
            else
               RAP = 1000000.0*sign(1.,PZ)
            endif
            IY = 1 + int(ABS(RAP) / DY)
            IF (IY.lt.1 .or.IY .GT. 50) GOTO 100
            dyp1(IY) = dyp1(IY) + 1.0
            DEYP1(IY) = DEYP1(IY) + XMT
            IF (ITYP .EQ. 21) THEN
               dyg1(IY) = dyg1(IY) + 1.0
               DEYG1(IY) = DEYG1(IY) + XMT
            END IF
 100        CONTINUE
            IMT = 1 + int(DXMT / DMT)
            IF (RAP .GT. YMAX .OR. RAP .LE. YMIN) GOTO 200
            IF (IMT .GT. 200) GOTO 200
            DMYP1(IMT) = DMYP1(IMT) + 1.0 / XMT
            IF (ITYP .EQ. 21) THEN
               DMYG1(IMT) = DMYG1(IMT) + 1.0 / XMT
            END IF
 200        CONTINUE
 1005    CONTINUE
 1006 CONTINUE
      DO 1008 I = 1, IHNT2(3)
         DO 1007 J = 1, NTJ(I)
            ITYP = KFTJ(I, J)
            PX = PJTX(I, J)
            PY = PJTY(I, J)
            PZ = PJTZ(I, J)
            PE = PJTE(I, J)
            PM = PJTM(I, J)
            XMT = SQRT(PX ** 2 + PY ** 2 + PM ** 2)
            DXMT = XMT - PM
            if(XMT.gt.0.) then
               RAP=asinh(PZ/XMT)
            else
               PRINT *, ' IN HJANA1 mt=0'
               RAP = 1000000.0*sign(1.,PZ)
            endif
            IY = 1 + int(ABS(RAP) / DY)
            IF (IY.lt.1 .or.IY .GT. 50) GOTO 300
            dyp1(IY) = dyp1(IY) + 1.0
            DEYP1(IY) = DEYP1(IY) + XMT
            IF (ITYP .EQ. 21) THEN
               dyg1(IY) = dyg1(IY) + 1.0
               DEYG1(IY) = DEYG1(IY) + XMT
            END IF
 300        CONTINUE
            IF (RAP .GT. YMAX .OR. RAP .LE. YMIN) GOTO 400
            IMT = 1 + int(DXMT / DMT)
            IF (IMT .GT. 200) GOTO 400
            DMYP1(IMT) = DMYP1(IMT) + 1.0 / XMT
            IF (ITYP .EQ. 21) THEN
               DMYG1(IMT) = DMYG1(IMT) + 1.0 / XMT
            END IF
 400        CONTINUE
 1007    CONTINUE
 1008 CONTINUE
      DO 1010 I = 1, NSG
         DO 1009 J = 1, NJSG(I)
            ITYP = K2SG(I, J)
            PX = PXSG(I, J)
            PY = PYSG(I, J)
            PZ = PZSG(I, J)
            PE = PESG(I, J)
            PM = PMSG(I, J)
            XMT = SQRT(PX ** 2 + PY ** 2 + PM ** 2)
            DXMT = XMT - PM
            if(XMT.gt.0.) then
               RAP=asinh(PZ/XMT)
            else
               PRINT *, ' IN HJANA1 mt=0'
               RAP = 1000000.0*sign(1.,PZ)
            endif
            IY = 1 + int(ABS(RAP) / DY)
            IF (IY.lt.1 .or.IY .GT. 50) GOTO 500
            dyp1(IY) = dyp1(IY) + 1.0
            DEYP1(IY) = DEYP1(IY) + XMT
            IF (ITYP .EQ. 21) THEN
               dyg1(IY) = dyg1(IY) + 1.0
               DEYG1(IY) = DEYG1(IY) + XMT
            END IF
 500        CONTINUE
            IF (RAP .GT. YMAX .OR. RAP .LE. YMIN) GOTO 600
            IMT = 1 + int(DXMT / DMT)
            IF (IMT .GT. 200) GOTO 600
            DMYP1(IMT) = DMYP1(IMT) + 1.0 / XMT
            IF (ITYP .EQ. 21) THEN
               DMYG1(IMT) = DMYG1(IMT) + 1.0 / XMT
            END IF
 600        CONTINUE
 1009    CONTINUE
 1010 CONTINUE
      DO 1011 I = 1, IHNT2(1)
         YR = SQRT(YP(1, I) ** 2 + YP(2, I) ** 2)
         IR = 1 + int(YR / DR)
         IF (IR .GT. 50 .or. IR .LT. 1) GOTO 601
         dnrpj1(IR) = dnrpj1(IR) + 1.0
         dnrtt1(IR) = dnrtt1(IR) + 1.0
 601     CONTINUE
 1011 CONTINUE
      DO 1012 I = 1, IHNT2(3)
         YR = SQRT(YT(1, I) ** 2 + YT(2, I) ** 2)
         IR = 1 + int(YR / DR)
         IF (IR .GT. 50 .or. IR .LT. 1) GOTO 602
         dnrtg1(IR) = dnrtg1(IR) + 1.0
         dnrtt1(IR) = dnrtt1(IR) + 1.0
 602     CONTINUE
 1012 CONTINUE
      DO 1013 I = 1, NSG
         Y1 = 0.5 * (YP(1, IASG(I, 1)) + YT(1, IASG(I, 2)))
         Y2 = 0.5 * (YP(2, IASG(I, 1)) + YT(2, IASG(I, 2)))
         YR = SQRT(Y1 ** 2 + Y2 ** 2)
         IR = 1 + int(YR / DR)
         IF (IR .GT. 50 .or. IR .LT. 1) GOTO 603
         dnrin1(IR) = dnrin1(IR) + 1.0
         dnrtt1(IR) = dnrtt1(IR) + 1.0
 603     CONTINUE
 1013 CONTINUE
      DO 1014 I = 1, MUL
         ITYP = ITYP0(I)
         PX = sngl(PX0(I))
         PY = sngl(PY0(I))
         PZ = sngl(PZ0(I))
         PE = sngl(E0(I))
         PM = sngl(XMASS0(I))
         XMT = SQRT(PX ** 2 + PY ** 2 + PM ** 2)
         DXMT = XMT - PM
         if(XMT.gt.0.) then
            RAP=asinh(PZ/XMT)
         else
            PRINT *, ' IN HJANA1 mt=0'
            RAP = 1000000.0*sign(1.,PZ)
         endif
         IY = 1 + int(ABS(RAP) / DY)
         IF (IY.lt.1 .or.IY .GT. 50) GOTO 700
         dyg1c(IY) = dyg1c(IY) + 1.0
         deyg1c(IY) = deyg1c(IY) + XMT
 700     CONTINUE
         IF (RAP .GT. YMAX .OR. RAP .LE. YMIN) GOTO 800
         IMT = 1 + int(DXMT / DMT)
         IF (IMT .GT. 50) GOTO 800
         dmyg1c(IMT) = dmyg1c(IMT) + 1.0 / XMT
 800     CONTINUE
 1014 CONTINUE
      DO 1016 I = 1, IHNT2(1)
         DO 1015 J = 1, NPJ(I)
            nsubp = nsubp + 1
            IF (KFPJ(I, J) .EQ. 21) nsubg = nsubg + 1
 1015    CONTINUE
 1016 CONTINUE
      DO 1018 I = 1, IHNT2(3)
         DO 1017 J = 1, NTJ(I)
            nsubp = nsubp + 1
            IF (KFTJ(I, J) .EQ. 21) nsubg = nsubg + 1
 1017    CONTINUE
 1018 CONTINUE
      DO 1020 I = 1, NSG
         DO 1019 J = 1, NJSG(I)
            nsubp = nsubp + 1
            IF (K2SG(I, J) .EQ. 21) nsubg = nsubg + 1
 1019    CONTINUE
 1020 CONTINUE
      NISG = NISG + NSG
      IF (IOUT .EQ. 1) THEN
      PRINT *, ' in HJANA1 '
      PRINT *, ' total number of partons = ', nsubp / IW
      PRINT *, ' total number of gluons = ', nsubg / IW
      PRINT *, ' number of independent strings = ', NISG / IW
      END IF
      RETURN
      END
