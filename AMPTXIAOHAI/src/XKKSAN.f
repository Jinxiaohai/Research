      SUBROUTINE XKKSAN(i1,i2,SRT,SIGKS1,SIGKS2,SIGKS3,SIGKS4,SIGK,prkk)
          PARAMETER (AKA=0.498, PIMASS=0.140, RHOM = 0.770,aks=0.895,
     & OMEGAM = 0.7819, ETAM = 0.5473)
      PARAMETER (MAXSTR=150001)
      COMMON  /CC/      E(MAXSTR)
      SAVE   
        S = SRT ** 2
       SIGKS1 = 1.E-08
       SIGKS2 = 1.E-08
       SIGKS3 = 1.E-08
       SIGKS4 = 1.E-08
        XPION0 = prkk
        XPION0 = XPION0/2
        PI2 = (S - (e(i1) + e(i2)) ** 2) * (S - (e(i1) - e(i2)) ** 2)
        SIGK = 1.E-08
        if(PI2 .le. 0.0) return
        XM1 = PIMASS
        XM2 = RHOM
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PI2 .GT. 0.0 .AND. PF2 .GT. 0.0) THEN
           SIGKS1 = 27.0 / 4.0 * PF2 / PI2 * XPION0
        END IF
        XM1 = PIMASS
        XM2 = OMEGAM
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PI2 .GT. 0.0 .AND. PF2 .GT. 0.0) THEN
           SIGKS2 = 9.0 / 4.0 * PF2 / PI2 * XPION0
        END IF
        XM1 = RHOM
        XM2 = ETAM
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           SIGKS3 = 9.0 / 4.0 * PF2 / PI2 * XPION0
        END IF
        XM1 = OMEGAM
        XM2 = ETAM
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           SIGKS4 = 3.0 / 4.0 * PF2 / PI2 * XPION0
        END IF
        SIGK=SIGKS1+SIGKS2+SIGKS3+SIGKS4
       RETURN
        END
