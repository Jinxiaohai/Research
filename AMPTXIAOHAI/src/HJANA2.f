      SUBROUTINE HJANA2
      PARAMETER (YMAX = 1.0, YMIN = -1.0)
      PARAMETER (DMT = 0.05, DY = 0.2)
      PARAMETER (DR = 0.2, DT = 0.2)
      PARAMETER (MAXPTN=400001)
      PARAMETER (MAXSTR=150001)
      DOUBLE PRECISION  PXSGS,PYSGS,PZSGS,PESGS,PMSGS,
     1     GXSGS,GYSGS,GZSGS,FTSGS
      DIMENSION dyp2(50), DMYP2(200), DEYP2(50)
      DIMENSION dyg2(50), DMYG2(200), DEYG2(50)
      DIMENSION SNYP2(50), SMYP2(200), SEYP2(50)
      DIMENSION SNYG2(50), SMYG2(200), SEYG2(50)
      DIMENSION dnrpj2(50), dnrtg2(50), dnrin2(50),
     &   dnrtt2(50)
      DIMENSION dtpj2(50), dttg2(50), dtin2(50),
     &   dttot2(50)
      DIMENSION dyg2c(50), dmyg2c(50), deyg2c(50)
      DIMENSION snrpj2(50), snrtg2(50), snrin2(50),
     &   snrtt2(50)
      DIMENSION stpj2(50), sttg2(50), stin2(50),
     &   sttot2(50)
      DIMENSION snyg2c(50), smyg2c(50), seyg2c(50)
      DOUBLE PRECISION  ATAUI, ZT1, ZT2, ZT3
      DOUBLE PRECISION  GX5, GY5, GZ5, FT5, PX5, PY5, PZ5, E5, XMASS5
      COMMON /PARA1/ MUL
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      COMMON /SREC2/ATAUI(MAXSTR),ZT1(MAXSTR),ZT2(MAXSTR),ZT3(MAXSTR)
      COMMON/HJJET1/NPJ(300),KFPJ(300,500),PJPX(300,500),
     &   PJPY(300,500),PJPZ(300,500),PJPE(300,500),
     &   PJPM(300,500),NTJ(300),KFTJ(300,500),
     &   PJTX(300,500),PJTY(300,500),PJTZ(300,500),
     &   PJTE(300,500),PJTM(300,500)
      COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &   K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &   PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
      COMMON /prec2/GX5(MAXPTN),GY5(MAXPTN),GZ5(MAXPTN),FT5(MAXPTN),
     &   PX5(MAXPTN), PY5(MAXPTN), PZ5(MAXPTN), E5(MAXPTN),
     &   XMASS5(MAXPTN), ITYP5(MAXPTN)
      COMMON /AREVT/ IAEVT, IARUN, MISS
      COMMON /AROUT/ IOUT
      common/anim/nevent,isoft,isflag,izpc
      COMMON/SOFT/PXSGS(MAXSTR,3),PYSGS(MAXSTR,3),PZSGS(MAXSTR,3),
     &     PESGS(MAXSTR,3),PMSGS(MAXSTR,3),GXSGS(MAXSTR,3),
     &     GYSGS(MAXSTR,3),GZSGS(MAXSTR,3),FTSGS(MAXSTR,3),
     &     K1SGS(MAXSTR,3),K2SGS(MAXSTR,3),NJSGS(MAXSTR)
      SAVE   
      DATA IW/0/
      IF (isevt .EQ. IAEVT .AND. isrun .EQ. IARUN) THEN
         DO 1001 I = 1, 200
            DMYP2(I) = SMYP2(I)
            DMYG2(I) = SMYG2(I)
 1001    CONTINUE
         DO 1002 I = 1, 50
            dyp2(I) = SNYP2(I)
            DEYP2(I) = SEYP2(I)
            dyg2(I) = SNYG2(I)
            DEYG2(I) = SEYG2(I)
            dnrpj2(I) = snrpj2(I)
            dnrtg2(I) = snrtg2(I)
            dnrin2(I) = snrin2(I)
            dnrtt2(I) = snrtt2(I)
            dtpj2(I) = stpj2(I)
            dttg2(I) = sttg2(I)
            dtin2(I) = stin2(I)
            dttot2(I) = sttot2(I)
            dyg2c(I) = snyg2c(I)
            dmyg2c(I) = smyg2c(I)
            deyg2c(I) = seyg2c(I)
 1002    CONTINUE
         nsubp = nsubpS
         nsubg = nsubgS
         NISG = NISGS
      ELSE
         DO 1003 I = 1, 200
            SMYP2(I) = DMYP2(I)
            SMYG2(I) = DMYG2(I)
 1003    CONTINUE
         DO 1004 I = 1, 50
            SNYP2(I) = dyp2(I)
            SEYP2(I) = DEYP2(I)
            SNYG2(I) = dyg2(I)
            SEYG2(I) = DEYG2(I)
            snrpj2(I) = dnrpj2(I)
            snrtg2(I) = dnrtg2(I)
            snrin2(I) = dnrin2(I)
            snrtt2(I) = dnrtt2(I)
            stpj2(I) = dtpj2(I)
            sttg2(I) = dttg2(I)
            stin2(I) = dtin2(I)
            sttot2(I) = dttot2(I)
            snyg2c(I) = dyg2c(I)
            smyg2c(I) = dmyg2c(I)
            seyg2c(I) = deyg2c(I)
 1004    CONTINUE
         nsubpS = nsubp
         nsubgS = nsubg
         NISGS = NISG
         isevt = IAEVT
         isrun = IARUN
         IW = IW + 1
      END IF
      if(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) goto 510
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
               PRINT *, ' IN HJANA2 mt=0'
               RAP = 1000000.0*sign(1.,PZ)
            endif
            IY = 1 + int(ABS(RAP) / DY)
            IF (IY.lt.1 .or.IY .GT. 50) GOTO 100
            dyp2(IY) = dyp2(IY) + 1.0
            DEYP2(IY) = DEYP2(IY) + XMT
            IF (ITYP .EQ. 21) THEN
               dyg2(IY) = dyg2(IY) + 1.0
               DEYG2(IY) = DEYG2(IY) + XMT
            END IF
 100        CONTINUE
            IF (RAP .GT. YMAX .OR. RAP .LE. YMIN) GOTO 200
            IMT = 1 + int(DXMT / DMT)
            IF (IMT .GT. 200) GOTO 200
            DMYP2(IMT) = DMYP2(IMT) + 1.0 / XMT
            IF (ITYP .EQ. 21) THEN
               DMYG2(IMT) = DMYG2(IMT) + 1.0 / XMT
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
               PRINT *, ' IN HJANA2 mt=0'
               RAP = 1000000.0*sign(1.,PZ)
            endif
            IY = 1 + int(ABS(RAP) / DY)
            IF (IY.lt.1 .or.IY .GT. 50) GOTO 300
            dyp2(IY) = dyp2(IY) + 1.0
            DEYP2(IY) = DEYP2(IY) + XMT
            IF (ITYP .EQ. 21) THEN
               dyg2(IY) = dyg2(IY) + 1.0
               DEYG2(IY) = DEYG2(IY) + XMT
            END IF
 300        CONTINUE
            IF (RAP .GT. YMAX .OR. RAP .LE. YMIN) GOTO 400
            IMT = 1 + int(DXMT / DMT)
            IF (IMT .GT. 200) GOTO 400
            DMYP2(IMT) = DMYP2(IMT) + 1.0 / XMT
            IF (ITYP .EQ. 21) THEN
               DMYG2(IMT) = DMYG2(IMT) + 1.0 / XMT
            END IF
 400        CONTINUE
 1007    CONTINUE
 1008 CONTINUE
 510  continue
      DO 1010 I = 1, NSG
         NJ=NJSG(I)
         if(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) NJ=NJSGS(I)
         DO 1009 J = 1, NJ
            ITYP = K2SG(I, J)
            PX = PXSG(I, J)
            PY = PYSG(I, J)
            PZ = PZSG(I, J)
            PE = PESG(I, J)
            PM = PMSG(I, J)
            if(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) then
               ITYP = K2SGS(I, J)
               PX = sngl(PXSGS(I, J))
               PY = sngl(PYSGS(I, J))
               PZ = sngl(PZSGS(I, J))
               PE = sngl(PESGS(I, J))
               PM = sngl(PMSGS(I, J))
            endif
            XMT = SQRT(PX ** 2 + PY ** 2 + PM ** 2)
            DXMT = XMT - PM
            if(XMT.gt.0.) then
               RAP=asinh(PZ/XMT)
            else
               PRINT *, ' IN HJANA2 mt=0'
               RAP = 1000000.0*sign(1.,PZ)
            endif
            IY = 1 + int(ABS(RAP) / DY)
            IF (IY.lt.1 .or.IY .GT. 50) GOTO 500
            dyp2(IY) = dyp2(IY) + 1.0
            DEYP2(IY) = DEYP2(IY) + XMT
            IF (ITYP .EQ. 21) THEN
               dyg2(IY) = dyg2(IY) + 1.0
               DEYG2(IY) = DEYG2(IY) + XMT
            END IF
 500        CONTINUE
            IF (RAP .GT. YMAX .OR. RAP .LE. YMIN) GOTO 600
            IMT = 1 + int(DXMT / DMT)
            IF (IMT .GT. 200) GOTO 600
            DMYP2(IMT) = DMYP2(IMT) + 1.0 / XMT
            IF (ITYP .EQ. 21) THEN
               DMYG2(IMT) = DMYG2(IMT) + 1.0 / XMT
            END IF
 600        CONTINUE
 1009    CONTINUE
 1010 CONTINUE
      if(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) goto 520
      DO 1011 I = 1, IHNT2(1)
         J = I
         YR = SQRT(sngl(ZT1(J) ** 2 + ZT2(J) ** 2))
         IR = 1 + int(YR / DR)
         IF (IR .GT. 50 .or. IR .LT. 1) GOTO 601
         dnrpj2(IR) = dnrpj2(IR) + 1.0
         dnrtt2(IR) = dnrtt2(IR) + 1.0
 601     CONTINUE
         IT = 1 + int(sngl(ATAUI(J)) / DT)
         IF (IT .GT. 50 .or. IT .LT. 1) GOTO 602
         dtpj2(IT) = dtpj2(IT) + 1.0
         dttot2(IT) = dttot2(IT) + 1.0
 602     CONTINUE
 1011 CONTINUE
      DO 1012 I = 1, IHNT2(3)
         J = I + IHNT2(1)
         YR = SQRT(sngl(ZT1(J) ** 2 + ZT2(J) ** 2))
         IR = 1 + int(YR / DR)
         IF (IR .GT. 50 .or. IR .LT. 1) GOTO 603
         dnrtg2(IR) = dnrtg2(IR) + 1.0
         dnrtt2(IR) = dnrtt2(IR) + 1.0
 603     CONTINUE
         IT = 1 + int(sngl(ATAUI(J)) / DT)
         IF (IT .GT. 50 .or. IT .LT. 1) GOTO 604
         dttg2(IT) = dttg2(IT) + 1.0
         dttot2(IT) = dttot2(IT) + 1.0
 604     CONTINUE
 1012 CONTINUE
 520  continue
      DO 1013 I = 1, NSG
         J = I + IHNT2(1) + IHNT2(3)
         if(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) J = I
         YR = SQRT(sngl(ZT1(J) ** 2 + ZT2(J) ** 2))
         IR = 1 + int(YR / DR)
         IF (IR .GT. 50 .or. IR .LT. 1) GOTO 605
         dnrin2(IR) = dnrin2(IR) + 1.0
         dnrtt2(IR) = dnrtt2(IR) + 1.0
 605     CONTINUE
         IT = 1 + int(sngl(ATAUI(J)) / DT)
         IF (IT .GT. 50 .or. IT .LT. 1) GOTO 606
         dtin2(IT) = dtin2(IT) + 1.0
         dttot2(IT) = dttot2(IT) + 1.0
 606     CONTINUE
 1013 CONTINUE
      DO 1014 I = 1, MUL
         ITYP = ITYP5(I)
         PX = sngl(PX5(I))
         PY = sngl(PY5(I))
         PZ = sngl(PZ5(I))
         PE = sngl(E5(I))
         PM = sngl(XMASS5(I))
         XMT = SQRT(PX ** 2 + PY ** 2 + PM ** 2)
         DXMT = XMT - PM
         if(XMT.gt.0.) then
            RAP=asinh(PZ/XMT)
         else
            PRINT *, ' IN HJANA2 mt=0'
            RAP = 1000000.0*sign(1.,PZ)
         endif
         IY = 1 + int(ABS(RAP) / DY)
         IF (IY.lt.1 .or.IY .GT. 50) GOTO 700
         dyg2c(IY) = dyg2c(IY) + 1.0
         deyg2c(IY) = deyg2c(IY) + XMT
 700     CONTINUE
         IF (RAP .GT. YMAX .OR. RAP .LE. YMIN) GOTO 800
         IMT = 1 + int(DXMT / DMT)
         IF (IMT .GT. 50) GOTO 800
         dmyg2c(IMT) = dmyg2c(IMT) + 1.0 / XMT
 800     CONTINUE
 1014 CONTINUE
      if(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) goto 530
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
 530  continue
      DO 1020 I = 1, NSG
         NJ=NJSG(I)
         if(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) NJ=NJSGS(I)
         DO 1019 J = 1, NJ
            nsubp = nsubp + 1
            if(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) then
               IF(K2SGS(I, J) .EQ. 21) nsubg = nsubg + 1
            else
               IF (K2SG(I, J) .EQ. 21) nsubg = nsubg + 1
            endif
 1019    CONTINUE
 1020 CONTINUE
      NISG = NISG + NSG
      IF (IOUT .EQ. 1) THEN
      PRINT *, ' in HJANA2 '
      PRINT *, ' total number of partons = ', nsubp / IW
      PRINT *, ' total number of gluons = ', nsubg / IW
      PRINT *, ' number of independent strings = ', NISG / IW
      END IF
      CALL HJAN2A
      CALL HJAN2B
      RETURN
      END
