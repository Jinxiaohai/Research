      SUBROUTINE XKKANN(SRT, XSK1, XSK2, XSK3, XSK4, XSK5,
     &     XSK6, XSK7, XSK8, XSK9, XSK10, XSK11, SIGK, rrkk)
      PARAMETER  (MAXSTR=150001, MAXX=20,  MAXZ=24)
          PARAMETER (AKA=0.498, PIMASS=0.140, RHOM = 0.770, 
     &     OMEGAM = 0.7819, ETAM = 0.5473, APHI=1.02)
      COMMON  /AA/ R(3,MAXSTR)
      COMMON /BB/  P(3,MAXSTR)
      COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
      COMMON  /DD/      RHO(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHOP(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHON(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ)
      SAVE   
        S = SRT ** 2
       SIGK = 1.E-08
        XSK1 = 0.0
        XSK2 = 0.0
        XSK3 = 0.0
        XSK4 = 0.0
        XSK5 = 0.0
        XSK6 = 0.0
        XSK7 = 0.0
        XSK8 = 0.0
        XSK9 = 0.0
        XSK10 = 0.0
        XSK11 = 0.0
        XPION0 = PIPIK(SRT)
        XPION0 = 2.0 * XPION0
        PI2 = S * (S - 4.0 * AKA ** 2)
         if(PI2 .le. 0.0)return
        XM1 = PIMASS
        XM2 = PIMASS
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XSK1 = 9.0 / 4.0 * PF2 / PI2 * XPION0
        END IF
        XM1 = PIMASS
        XM2 = ETAM
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XSK4 = 3.0 / 4.0 * PF2 / PI2 * XPION0
        END IF
        XM1 = ETAM
        XM2 = ETAM
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XSK10 = 1.0 / 4.0 * PF2 / PI2 * XPION0
        END IF
        XPION0 = rrkk
        XM1 = RHOM
        XM2 = RHOM
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XSK5 = 81.0 / 4.0 * PF2 / PI2 * XPION0
        END IF
        XM1 = RHOM
        XM2 = OMEGAM
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XSK6 = 27.0 / 4.0 * PF2 / PI2 * XPION0
        END IF
        XM1 = OMEGAM
        XM2 = OMEGAM
        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
        IF (PF2 .GT. 0.0) THEN
           XSK8 = 9.0 / 4.0 * PF2 / PI2 * XPION0
        END IF
          fwdp = 1.68*(aphi**2-4.*aka**2)**1.5/6./aphi/aphi     
          scheck=srt**2-4.0*aka**2
          if(scheck.le.0) then
             write(99,*) 'scheck47: ', scheck
             stop
          endif
          pkaon=0.5*sqrt(scheck)
          XSK11 = 30.*3.14159*0.1973**2*(aphi*fwdp)**2/
     &             ((srt**2-aphi**2)**2+(aphi*fwdp)**2)/pkaon**2
        SIGK = XSK1 + XSK2 + XSK3 + XSK4 + XSK5 + 
     &     XSK6 + XSK7 + XSK8 + XSK9 + XSK10 + XSK11
       RETURN
        END
