        SUBROUTINE XKHYPE(I1, I2, SRT, XKY1, XKY2, XKY3, XKY4, XKY5,
     &     XKY6, XKY7, XKY8, XKY9, XKY10, XKY11, XKY12, XKY13,
     &     XKY14, XKY15, XKY16, XKY17, SIGK)
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,AMRHO=0.769,AMOMGA=0.782,APHI=1.02,
     2  AP2=0.13957,AM0=1.232,PI=3.1415926,CUTOFF=1.8966,AVMASS=0.9383)
          parameter (pimass=0.140, AMETA = 0.5473, aka=0.498,
     &     aml=1.116,ams=1.193, AM1440 = 1.44, AM1535 = 1.535)
        COMMON  /EE/ID(MAXSTR), LB(MAXSTR)
      SAVE   
        S = SRT ** 2
       SIGK=1.E-08 
        XKY1 = 0.0
        XKY2 = 0.0
        XKY3 = 0.0
        XKY4 = 0.0
        XKY5 = 0.0
        XKY6 = 0.0
        XKY7 = 0.0
        XKY8 = 0.0
        XKY9 = 0.0
        XKY10 = 0.0
        XKY11 = 0.0
        XKY12 = 0.0
        XKY13 = 0.0
        XKY14 = 0.0
        XKY15 = 0.0
        XKY16 = 0.0
        XKY17 = 0.0
        LB1 = LB(I1)
        LB2 = LB(I2)
        IF (iabs(LB1) .EQ. 14 .OR. iabs(LB2) .EQ. 14) THEN
           XKAON0 = PNLKA(SRT)
           XKAON0 = 2.0 * XKAON0
           PI2 = (S - (AML + AKA) ** 2) * (S - (AML - AKA) ** 2)
        ELSE
           XKAON0 = PNSKA(SRT)
           XKAON0 = 2.0 * XKAON0
           PI2 = (S - (AMS + AKA) ** 2) * (S - (AMS - AKA) ** 2)
        END IF
          if(PI2 .le. 0.0)return
        XM1 = PIMASS
        XM2 = AMP
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XKY1 = 3.0 * PF2 / PI2 * XKAON0
        END IF
        XM1 = PIMASS
        XM2 = AM0
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XKY2 = 12.0 * PF2 / PI2 * XKAON0
        END IF
        XM1 = PIMASS
        XM2 = AM1440
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XKY3 = 3.0 * PF2 / PI2 * XKAON0
        END IF
        XM1 = PIMASS
        XM2 = AM1535
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XKY4 = 3.0 * PF2 / PI2 * XKAON0
        END IF
        XM1 = AMRHO
        XM2 = AMP
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XKY5 = 9.0 * PF2 / PI2 * XKAON0
        END IF
        XM1 = AMRHO
        XM2 = AM0
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XKY6 = 36.0 * PF2 / PI2 * XKAON0
        END IF
        XM1 = AMRHO
        XM2 = AM1440
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XKY7 = 9.0 * PF2 / PI2 * XKAON0
        END IF
        XM1 = AMRHO
        XM2 = AM1535
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XKY8 = 9.0 * PF2 / PI2 * XKAON0
        END IF
        XM1 = AMOMGA
        XM2 = AMP
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XKY9 = 3.0 * PF2 / PI2 * XKAON0
        END IF
        XM1 = AMOMGA
        XM2 = AM0
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XKY10 = 12.0 * PF2 / PI2 * XKAON0
        END IF
        XM1 = AMOMGA
        XM2 = AM1440
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XKY11 = 3.0 * PF2 / PI2 * XKAON0
        END IF
        XM1 = AMOMGA
        XM2 = AM1535
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XKY12 = 3.0 * PF2 / PI2 * XKAON0
        END IF
        XM1 = AMETA
        XM2 = AMP
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XKY13 = 1.0 * PF2 / PI2 * XKAON0
        END IF
        XM1 = AMETA
        XM2 = AM0
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XKY14 = 4.0 * PF2 / PI2 * XKAON0
        END IF
        XM1 = AMETA
        XM2 = AM1440
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XKY15 = 1.0 * PF2 / PI2 * XKAON0
        END IF
        XM1 = AMETA
        XM2 = AM1535
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XKY16 = 1.0 * PF2 / PI2 * XKAON0
        END IF
        if(lb1.eq.14 .or. lb2.eq.14)then
         if(srt .gt. (aphi+amn))then
           srrt = srt - (aphi+amn)
           sig = 1.715/((srrt+3.508)**2-12.138)
         XM1 = AMN
         XM2 = APHI
         PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
         XKY17 = 3.0 * PF2 / PI2 * SIG/10.
        endif
       endif
       IF ((iabs(LB1) .GE. 15 .AND. iabs(LB1) .LE. 17) .OR. 
     &     (iabs(LB2) .GE. 15 .AND. iabs(LB2) .LE. 17)) THEN
           DDF = 3.0
           XKY1 = XKY1 / DDF
           XKY2 = XKY2 / DDF
           XKY3 = XKY3 / DDF
           XKY4 = XKY4 / DDF
           XKY5 = XKY5 / DDF
           XKY6 = XKY6 / DDF
           XKY7 = XKY7 / DDF
           XKY8 = XKY8 / DDF
           XKY9 = XKY9 / DDF
           XKY10 = XKY10/ DDF
           XKY11 = XKY11 / DDF
           XKY12 = XKY12 / DDF
           XKY13 = XKY13 / DDF
           XKY14 = XKY14 / DDF
           XKY15 = XKY15 / DDF
           XKY16 = XKY16 / DDF
        END IF
        SIGK = XKY1 + XKY2 + XKY3 + XKY4 +
     &       XKY5 + XKY6 + XKY7 + XKY8 +
     &       XKY9 + XKY10 + XKY11 + XKY12 +
     &       XKY13 + XKY14 + XKY15 + XKY16 + XKY17
       RETURN
       END
