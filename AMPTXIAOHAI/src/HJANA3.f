      SUBROUTINE HJANA3
      PARAMETER (MAXSTR=150001, MAXR=1)
      PARAMETER (YMIN = -1.0, YMAX = 1.0)
      PARAMETER (DMT = 0.05, DY = 0.2)
      DOUBLE PRECISION v2i,eti,xmulti,v2mi,s2mi,xmmult,
     1     v2bi,s2bi,xbmult
      DIMENSION dndyh3(50), DMYH3(50), DEYH3(50)
      COMMON /RUN/ NUM
      COMMON /ARERC1/MULTI1(MAXR)
      COMMON /ARPRC1/ITYP1(MAXSTR, MAXR),
     &     GX1(MAXSTR, MAXR), GY1(MAXSTR, MAXR), GZ1(MAXSTR, MAXR), 
     &     FT1(MAXSTR, MAXR),
     &     PX1(MAXSTR, MAXR), PY1(MAXSTR, MAXR), PZ1(MAXSTR, MAXR),
     &     EE1(MAXSTR, MAXR), XM1(MAXSTR, MAXR)
      COMMON /AROUT/ IOUT
      COMMON/iflow/v2i,eti,xmulti,v2mi,s2mi,xmmult,v2bi,s2bi,xbmult
      SAVE   
      DATA IW/0/
      IW = IW + 1
      DO 1002 J = 1, NUM
         DO 1001 I = 1, MULTI1(J)
            ITYP = ITYP1(I, J)
            IF (ITYP .GT. -100 .AND. ITYP .LT. 100) GOTO 200
            PX = PX1(I, J)
            PY = PY1(I, J)
            PZ = PZ1(I, J)
            EE = EE1(I, J)
            XM = XM1(I, J)
            XMT = SQRT(PX ** 2 + PY ** 2 + XM ** 2)
            DXMT = XMT - XM
            if(XMT.gt.0.) then
               Y=asinh(PZ/XMT)
            else
               PRINT *, ' IN HJANA3 mt=0'
               Y = 1000000.0*sign(1.,PZ)
            endif
            IY = 1 + int(Y/DY)
            IF (IY.lt.1 .or.IY .GT. 50) GOTO 100
            dndyh3(IY) = dndyh3(IY) + 1.0
            DEYH3(IY) = DEYH3(IY) + XMT
 100        CONTINUE
            IF (Y. LT. YMIN .OR. Y .GE. YMAX) GOTO 200
            IMT = 1 + int(DXMT / DMT)
            IF (IMT .GT. 50) GOTO 200
            DMYH3(IMT) = DMYH3(IMT) + 1.0 / XMT
 200        CONTINUE
 1001    CONTINUE
 1002 CONTINUE
      RETURN
      END
