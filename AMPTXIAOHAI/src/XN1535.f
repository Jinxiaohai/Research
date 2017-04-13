        REAL FUNCTION XN1535(I1,I2,LA)
        PARAMETER (MAXSTR=150001,MAXR=1,
     1  AMN=0.939457,AMP=0.93828,ETAM=0.5475,
     2  AP1=0.13496,AP2=0.13957,AM0=1.232,PI=3.1415926)
        double precision e10,e20,scheck,p1,p2,p3
        COMMON /AA/ R(3,MAXSTR)
        COMMON /BB/ P(3,MAXSTR)
        COMMON /CC/ E(MAXSTR)
        COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
        COMMON   /RUN/NUM
        COMMON   /PA/RPION(3,MAXSTR,MAXR)
        COMMON   /PB/PPION(3,MAXSTR,MAXR)
        COMMON   /PC/EPION(MAXSTR,MAXR)
        COMMON   /PD/LPION(MAXSTR,MAXR)
      SAVE   
        AVMASS=0.5*(AMN+AMP)
        AVPI=(2.*AP2+AP1)/3.
        E10=dSQRT(dble(E(I1))**2+dble(P(1,I1))**2
     1     +dble(P(2,I1))**2+dble(P(3,I1))**2)
        E20=dSQRT(dble(E(I2))**2+dble(P(1,I2))**2
     1       +dble(P(2,I2))**2+dble(P(3,I2))**2)
        p1=dble(P(1,I1))+dble(P(1,I2))
        p2=dble(P(2,I1))+dble(P(2,I2))
        p3=dble(P(3,I1))+dble(P(3,I2))
        scheck=(E10+E20)**2-p1**2-p2**2-p3**2
        if(scheck.lt.0) then
           write(99,*) 'scheck21: ', scheck
           scheck=0.d0
        endif
        DM=SQRT(sngl(scheck))
        IF(DM.LE.1.1) THEN
        XN1535=1.E-06
        RETURN
        ENDIF
        GAM=W1535(DM)
       GAM0=0.15
        F1=0.25*GAM0**2/(0.25*GAM**2+(DM-1.535)**2)
        IF(LA.EQ.1)THEN
       XMAX=11.3
        ELSE
       XMAX=74.
        ENDIF
        XN1535=F1*XMAX/10.
        RETURN
        END
