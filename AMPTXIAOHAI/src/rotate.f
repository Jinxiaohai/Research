      SUBROUTINE rotate(PX0,PY0,PZ0,px,py,pz)
      SAVE   
      PR0  = SQRT( PX0**2 + PY0**2 + PZ0**2 )
      IF(PR0.EQ.0)PR0=0.00000001
      C2  = PZ0 / PR0
      IF(PX0 .EQ. 0.0 .AND. PY0 .EQ. 0.0) THEN
        T2 = 0.0
      ELSE
        T2=ATAN2(PY0,PX0)
      END IF
      scheck=1.0 - C2**2
      if(scheck.lt.0) then
         write(99,*) 'scheck45: ', scheck
         scheck=0.
      endif
      S2=sqrt(scheck)
      CT2  = COS(T2)
      ST2  = SIN(T2)
      PR=SQRT(PX**2+PY**2+PZ**2)
      IF(PR.EQ.0)PR=0.0000001
      C1=PZ/PR
      IF(PX.EQ.0.AND.PY.EQ.0)THEN
      T1=0.
      ELSE
      T1=ATAN2(PY,PX)
      ENDIF
      scheck=1.0 - C1**2
      if(scheck.lt.0) then
         write(99,*) 'scheck46: ', scheck
         scheck=0.
      endif
      S1=sqrt(scheck)
      CT1  = COS(T1)
      ST1  = SIN(T1)
      SS   = C2 * S1 * CT1  +  S2 * C1
      PX   = PR * ( SS*CT2 - S1*ST1*ST2 )
      PY   = PR * ( SS*ST2 + S1*ST1*CT2 )
      PZ   = PR * ( C1*C2 - S1*S2*CT1 )
      RETURN
      END
