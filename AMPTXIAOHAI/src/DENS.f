      SUBROUTINE DENS(IPOT,MASS,NUM,NESC)
      PARAMETER     (MAXSTR= 150001,MAXR=1)
      PARAMETER     (MAXX   =    20,  MAXZ  =    24)
      dimension pxl(-maxx:maxx,-maxx:maxx,-maxz:maxz),
     1          pyl(-maxx:maxx,-maxx:maxx,-maxz:maxz),
     2          pzl(-maxx:maxx,-maxx:maxx,-maxz:maxz)
      COMMON  /AA/      R(3,MAXSTR)
      COMMON  /BB/      P(3,MAXSTR)
      COMMON  /CC/      E(MAXSTR)
      COMMON  /DD/      RHO(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHOP(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHON(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ)
      COMMON  /DDpi/    piRHO(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ)
      COMMON  /EE/      ID(MAXSTR),LB(MAXSTR)
      common  /ss/  inout(20)
      COMMON  /RR/  MASSR(0:MAXR)
      common  /tt/  PEL(-maxx:maxx,-maxx:maxx,-maxz:maxz)
     &,rxy(-maxx:maxx,-maxx:maxx,-maxz:maxz)
      common  /bbb/ bxx(-maxx:maxx,-maxx:maxx,-maxz:maxz),
     &byy(-maxx:maxx,-maxx:maxx,-maxz:maxz),
     &bzz(-maxx:maxx,-maxx:maxx,-maxz:maxz)
      real zet(-45:45)
      SAVE   
      data zet /
     4     1.,0.,0.,0.,0.,
     3     1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     2     -1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     1     0.,0.,0.,-1.,0.,1.,0.,-1.,0.,-1.,
     s     0.,-2.,-1.,0.,1.,0.,0.,0.,0.,-1.,
     e     0.,
     s     1.,0.,-1.,0.,1.,-1.,0.,1.,2.,0.,
     1     1.,0.,1.,0.,-1.,0.,1.,0.,0.,0.,
     2     -1.,0.,1.,0.,-1.,0.,1.,0.,0.,1.,
     3     0.,0.,0.,0.,0.,0.,0.,0.,0.,-1.,
     4     0.,0.,0.,0.,-1./
      DO 300 IZ = -MAXZ,MAXZ
        DO 200 IY = -MAXX,MAXX
          DO 100 IX = -MAXX,MAXX
            RHO(IX,IY,IZ) = 0.0
            RHOn(IX,IY,IZ) = 0.0
            RHOp(IX,IY,IZ) = 0.0
            piRHO(IX,IY,IZ) = 0.0
           pxl(ix,iy,iz) = 0.0
           pyl(ix,iy,iz) = 0.0
           pzl(ix,iy,iz) = 0.0
           pel(ix,iy,iz) = 0.0
           bxx(ix,iy,iz) = 0.0
           byy(ix,iy,iz) = 0.0
           bzz(ix,iy,iz) = 0.0
  100     CONTINUE
  200   CONTINUE
  300 CONTINUE
      NESC  = 0
      BIG   = 1.0 / ( 3.0 * FLOAT(NUM) )
      SMALL = 1.0 / ( 9.0 * FLOAT(NUM) )
      MSUM=0
      DO 400 IRUN = 1,NUM
      MSUM=MSUM+MASSR(IRUN-1)
      DO 400 J=1,MASSr(irun)
      I=J+MSUM
        IX = NINT( R(1,I) )
        IY = NINT( R(2,I) )
        IZ = NINT( R(3,I) )
        IF( IX .LE. -MAXX .OR. IX .GE. MAXX .OR.
     &      IY .LE. -MAXX .OR. IY .GE. MAXX .OR.
     &      IZ .LE. -MAXZ .OR. IZ .GE. MAXZ )    THEN
          NESC = NESC + 1
        ELSE
          if(j.gt.mass)go to 30
          RHO(IX,  IY,  IZ  ) = RHO(IX,  IY,  IZ  ) + BIG
          RHO(IX+1,IY,  IZ  ) = RHO(IX+1,IY,  IZ  ) + SMALL
          RHO(IX-1,IY,  IZ  ) = RHO(IX-1,IY,  IZ  ) + SMALL
          RHO(IX,  IY+1,IZ  ) = RHO(IX,  IY+1,IZ  ) + SMALL
          RHO(IX,  IY-1,IZ  ) = RHO(IX,  IY-1,IZ  ) + SMALL
          RHO(IX,  IY,  IZ+1) = RHO(IX,  IY,  IZ+1) + SMALL
          RHO(IX,  IY,  IZ-1) = RHO(IX,  IY,  IZ-1) + SMALL
         IF(ZET(LB(I)).NE.0)THEN
          RHOP(IX,  IY,  IZ  ) = RHOP(IX,  IY,  IZ  ) + BIG
          RHOP(IX+1,IY,  IZ  ) = RHOP(IX+1,IY,  IZ  ) + SMALL
          RHOP(IX-1,IY,  IZ  ) = RHOP(IX-1,IY,  IZ  ) + SMALL
          RHOP(IX,  IY+1,IZ  ) = RHOP(IX,  IY+1,IZ  ) + SMALL
          RHOP(IX,  IY-1,IZ  ) = RHOP(IX,  IY-1,IZ  ) + SMALL
          RHOP(IX,  IY,  IZ+1) = RHOP(IX,  IY,  IZ+1) + SMALL
          RHOP(IX,  IY,  IZ-1) = RHOP(IX,  IY,  IZ-1) + SMALL
         go to 40
         ENDIF
         IF(ZET(LB(I)).EQ.0)THEN
          RHON(IX,  IY,  IZ  ) = RHON(IX,  IY,  IZ  ) + BIG
          RHON(IX+1,IY,  IZ  ) = RHON(IX+1,IY,  IZ  ) + SMALL
          RHON(IX-1,IY,  IZ  ) = RHON(IX-1,IY,  IZ  ) + SMALL
          RHON(IX,  IY+1,IZ  ) = RHON(IX,  IY+1,IZ  ) + SMALL
          RHON(IX,  IY-1,IZ  ) = RHON(IX,  IY-1,IZ  ) + SMALL
          RHON(IX,  IY,  IZ+1) = RHON(IX,  IY,  IZ+1) + SMALL
          RHON(IX,  IY,  IZ-1) = RHON(IX,  IY,  IZ-1) + SMALL
         go to 40
          END IF
30              piRHO(IX,  IY,  IZ  ) = piRHO(IX,  IY,  IZ  ) + BIG
          piRHO(IX+1,IY,  IZ  ) = piRHO(IX+1,IY,  IZ  ) + SMALL
          piRHO(IX-1,IY,  IZ  ) = piRHO(IX-1,IY,  IZ  ) + SMALL
          piRHO(IX,  IY+1,IZ  ) = piRHO(IX,  IY+1,IZ  ) + SMALL
          piRHO(IX,  IY-1,IZ  ) = piRHO(IX,  IY-1,IZ  ) + SMALL
          piRHO(IX,  IY,  IZ+1) = piRHO(IX,  IY,  IZ+1) + SMALL
          piRHO(IX,  IY,  IZ-1) = piRHO(IX,  IY,  IZ-1) + SMALL
40       pxl(ix,iy,iz)=pxl(ix,iy,iz)+p(1,I)*BIG
       pxl(ix+1,iy,iz)=pxl(ix+1,iy,iz)+p(1,I)*SMALL
       pxl(ix-1,iy,iz)=pxl(ix-1,iy,iz)+p(1,I)*SMALL
       pxl(ix,iy+1,iz)=pxl(ix,iy+1,iz)+p(1,I)*SMALL
       pxl(ix,iy-1,iz)=pxl(ix,iy-1,iz)+p(1,I)*SMALL
       pxl(ix,iy,iz+1)=pxl(ix,iy,iz+1)+p(1,I)*SMALL
       pxl(ix,iy,iz-1)=pxl(ix,iy,iz-1)+p(1,I)*SMALL
       pYl(ix,iy,iz)=pYl(ix,iy,iz)+p(2,I)*BIG
       pYl(ix+1,iy,iz)=pYl(ix+1,iy,iz)+p(2,I)*SMALL
       pYl(ix-1,iy,iz)=pYl(ix-1,iy,iz)+p(2,I)*SMALL
       pYl(ix,iy+1,iz)=pYl(ix,iy+1,iz)+p(2,I)*SMALL
       pYl(ix,iy-1,iz)=pYl(ix,iy-1,iz)+p(2,I)*SMALL
       pYl(ix,iy,iz+1)=pYl(ix,iy,iz+1)+p(2,I)*SMALL
       pYl(ix,iy,iz-1)=pYl(ix,iy,iz-1)+p(2,I)*SMALL
       pZl(ix,iy,iz)=pZl(ix,iy,iz)+p(3,I)*BIG
       pZl(ix+1,iy,iz)=pZl(ix+1,iy,iz)+p(3,I)*SMALL
       pZl(ix-1,iy,iz)=pZl(ix-1,iy,iz)+p(3,I)*SMALL
       pZl(ix,iy+1,iz)=pZl(ix,iy+1,iz)+p(3,I)*SMALL
       pZl(ix,iy-1,iz)=pZl(ix,iy-1,iz)+p(3,I)*SMALL
       pZl(ix,iy,iz+1)=pZl(ix,iy,iz+1)+p(3,I)*SMALL
       pZl(ix,iy,iz-1)=pZl(ix,iy,iz-1)+p(3,I)*SMALL
       pel(ix,iy,iz)=pel(ix,iy,iz)
     1     +sqrt(e(I)**2+p(1,i)**2+p(2,I)**2+p(3,I)**2)*BIG
       pel(ix+1,iy,iz)=pel(ix+1,iy,iz)
     1     +sqrt(e(I)**2+p(1,i)**2+p(2,I)**2+p(3,I)**2)*SMALL
       pel(ix-1,iy,iz)=pel(ix-1,iy,iz)
     1     +sqrt(e(I)**2+p(1,i)**2+p(2,I)**2+p(3,I)**2)*SMALL
       pel(ix,iy+1,iz)=pel(ix,iy+1,iz)
     1     +sqrt(e(I)**2+p(1,i)**2+p(2,I)**2+p(3,I)**2)*SMALL
       pel(ix,iy-1,iz)=pel(ix,iy-1,iz)
     1     +sqrt(e(I)**2+p(1,i)**2+p(2,I)**2+p(3,I)**2)*SMALL
       pel(ix,iy,iz+1)=pel(ix,iy,iz+1)
     1     +sqrt(e(I)**2+p(1,i)**2+p(2,I)**2+p(3,I)**2)*SMALL
       pel(ix,iy,iz-1)=pel(ix,iy,iz-1)
     1     +sqrt(e(I)**2+p(1,i)**2+p(2,I)**2+p(3,I)**2)*SMALL
        END IF
  400 CONTINUE
      DO 301 IZ = -MAXZ,MAXZ
        DO 201 IY = -MAXX,MAXX
          DO 101 IX = -MAXX,MAXX
      IF((RHO(IX,IY,IZ).EQ.0).OR.(PEL(IX,IY,IZ).EQ.0))
     1GO TO 101
      SMASS2=PEL(IX,IY,IZ)**2-PXL(IX,IY,IZ)**2
     1-PYL(IX,IY,IZ)**2-PZL(IX,IY,IZ)**2
       IF(SMASS2.LE.0)SMASS2=1.E-06
       SMASS=SQRT(SMASS2)
           IF(SMASS.EQ.0.)SMASS=1.e-06
           GAMMA=PEL(IX,IY,IZ)/SMASS
           if(gamma.eq.0)go to 101
       bxx(ix,iy,iz)=pxl(ix,iy,iz)/pel(ix,iy,iz)                  
       byy(ix,iy,iz)=pyl(ix,iy,iz)/pel(ix,iy,iz)       
       bzz(ix,iy,iz)=pzl(ix,iy,iz)/pel(ix,iy,iz)                  
            RHO(IX,IY,IZ) = RHO(IX,IY,IZ)/GAMMA
            RHOn(IX,IY,IZ) = RHOn(IX,IY,IZ)/GAMMA
            RHOp(IX,IY,IZ) = RHOp(IX,IY,IZ)/GAMMA
            piRHO(IX,IY,IZ) = piRHO(IX,IY,IZ)/GAMMA
            pEL(IX,IY,IZ) = pEL(IX,IY,IZ)/(GAMMA**2)
           rho0=0.163
           IF(IPOT.EQ.0)THEN
           U=0
           GO TO 70
           ENDIF
           IF(IPOT.EQ.1.or.ipot.eq.6)THEN
           A=-0.1236
           B=0.0704
           S=2
           GO TO 60
           ENDIF
           IF(IPOT.EQ.2.or.ipot.eq.7)THEN
           A=-0.218
           B=0.164
           S=4./3.
           ENDIF
           IF(IPOT.EQ.3)THEN
           a=-0.3581
           b=0.3048
           S=1.167
           GO TO 60
           ENDIF
           IF(IPOT.EQ.4)THEN
           denr=rho(ix,iy,iz)/rho0         
           b=0.3048
           S=1.167
           if(denr.le.4.or.denr.gt.7)then
           a=-0.3581
           else
           a=-b*denr**(1./6.)-2.*0.036/3.*denr**(-0.333)
           endif
           GO TO 60
           ENDIF
60           U = 0.5*A*RHO(IX,IY,IZ)**2/RHO0 
     1        + B/(1+S) * (RHO(IX,IY,IZ)/RHO0)**S*RHO(IX,IY,IZ)  
70           PEL(IX,IY,IZ)=PEL(IX,IY,IZ)+U
  101     CONTINUE
  201   CONTINUE
  301 CONTINUE
      RETURN
      END
