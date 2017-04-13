      SUBROUTINE INIT(MINNUM,MAXNUM,NUM,RADIUS,X0,Z0,P0,
     &                GAMMA,ISEED,MASS,IOPT)
      PARAMETER     (MAXSTR=150001,  AMU   = 0.9383)
      PARAMETER     (MAXX   =   20,  MAXZ  =    24)
      PARAMETER     (PI=3.1415926)
      REAL              PTOT(3)
      COMMON  /AA/      R(3,MAXSTR)
      COMMON  /BB/      P(3,MAXSTR)
      COMMON  /CC/      E(MAXSTR)
      COMMON  /DD/      RHO(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHOP(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHON(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ)
      COMMON  /EE/      ID(MAXSTR),LB(MAXSTR)
      common  /ss/      inout(20)
      COMMON/RNDF77/NSEED
      SAVE   
      IF (P0 .NE. 0.) THEN
        SIGN = P0 / ABS(P0)
      ELSE
        SIGN = 0.
      END IF
      scheck=GAMMA**2-1.
      if(scheck.lt.0) then
         write(99,*) 'scheck10: ', scheck
         scheck=0.
      endif
      BETA=SIGN*SQRT(scheck)/GAMMA
      IF (MINNUM .EQ. 1) THEN
        IDNUM = 1
      ELSE
        IDNUM = -1
      END IF
      DO 400 IRUN = 1,NUM
        DO 100 I = MINNUM+(IRUN-1)*MASS,MAXNUM+(IRUN-1)*MASS
          ID(I) = IDNUM
          E(I)  = AMU
  100   CONTINUE
        DO 300 I = MINNUM+(IRUN-1)*MASS,MAXNUM+(IRUN-1)*MASS
  200     CONTINUE
            X = 1.0 - 2.0 * RANART(NSEED)
            Y = 1.0 - 2.0 * RANART(NSEED)
            Z = 1.0 - 2.0 * RANART(NSEED)
          IF ((X*X+Y*Y+Z*Z) .GT. 1.0) GOTO 200
          R(1,I) = X * RADIUS
          R(2,I) = Y * RADIUS
          R(3,I) = Z * RADIUS
  300   CONTINUE
  400 CONTINUE
      IF (IOPT .NE. 3) THEN
        RHOW0  = 0.168
        DO 1000 IRUN = 1,NUM
          DO 600 I = MINNUM+(IRUN-1)*MASS,MAXNUM+(IRUN-1)*MASS
  500       CONTINUE
              PX = 1.0 - 2.0 * RANART(NSEED)
              PY = 1.0 - 2.0 * RANART(NSEED)
              PZ = 1.0 - 2.0 * RANART(NSEED)
            IF (PX*PX+PY*PY+PZ*PZ .GT. 1.0) GOTO 500
            RDIST  = SQRT( R(1,I)**2 + R(2,I)**2 + R(3,I)**2 )
            RHOWS  = RHOW0 / (  1.0 + EXP( (RDIST-RADIUS) / 0.55 )  )
            PFERMI = 0.197 * (1.5 * PI*PI * RHOWS)**(1./3.)
            IF(IOPT.EQ.2) PFERMI=0.27
           if(iopt.eq.4) pfermi=0.
            P(1,I) = PFERMI * PX
            P(2,I) = PFERMI * PY
            P(3,I) = PFERMI * PZ
  600     CONTINUE
          DO 700 IDIR = 1,3
            PTOT(IDIR) = 0.0
  700     CONTINUE
          NPART = 0
          DO 900 I = MINNUM+(IRUN-1)*MASS,MAXNUM+(IRUN-1)*MASS
            NPART = NPART + 1
            DO 800 IDIR = 1,3
              PTOT(IDIR) = PTOT(IDIR) + P(IDIR,I)
  800       CONTINUE
  900     CONTINUE
          DO 950 I = MINNUM+(IRUN-1)*MASS,MAXNUM+(IRUN-1)*MASS
            DO 925 IDIR = 1,3
              P(IDIR,I) = P(IDIR,I) - PTOT(IDIR) / FLOAT(NPART)
  925       CONTINUE
            IF ((IOPT .EQ. 1).or.(iopt.eq.2)) THEN
              EPART = SQRT(P(1,I)**2+P(2,I)**2+P(3,I)**2+AMU**2)
              P(3,I) = GAMMA*(P(3,I) + BETA*EPART)
            ELSE
              P(3,I) = P(3,I) + P0
            END IF
  950     CONTINUE
 1000   CONTINUE
      ELSE
        DO 1200 IRUN = 1,NUM
          DO 1100 I = MINNUM+(IRUN-1)*MASS,MAXNUM+(IRUN-1)*MASS
            P(1,I) = 0.0
            P(2,I) = 0.0
            P(3,I) = P0
 1100     CONTINUE
 1200   CONTINUE
      END IF
      DO 1400 IRUN = 1,NUM
        DO 1300 I = MINNUM+(IRUN-1)*MASS,MAXNUM+(IRUN-1)*MASS
          R(1,I) = R(1,I) + X0
          R(3,I) = (R(3,I)+Z0)/ GAMMA 
 1300   CONTINUE
 1400 CONTINUE
      RETURN
      END
