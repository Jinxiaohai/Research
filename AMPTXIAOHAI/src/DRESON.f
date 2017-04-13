        SUBROUTINE DRESON(I1,I2)
        PARAMETER (MAXSTR=150001,MAXR=1,
     1  AMN=0.939457,AMP=0.93828,
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
        E10=dSQRT(dble(E(I1))**2+dble(P(1,I1))**2
     1     +dble(P(2,I1))**2+dble(P(3,I1))**2)
        E20=dSQRT(dble(E(I2))**2+dble(P(1,I2))**2
     1       +dble(P(2,I2))**2+dble(P(3,I2))**2)
        p1=dble(P(1,I1))+dble(P(1,I2))
        p2=dble(P(2,I1))+dble(P(2,I2))
        p3=dble(P(3,I1))+dble(P(3,I2))
        IF(iabs(LB(I2)) .EQ. 1 .OR. iabs(LB(I2)) .EQ. 2 .OR.
     &     (iabs(LB(I2)) .GE. 6 .AND. iabs(LB(I2)) .LE. 17)) THEN
        E(I1)=0.
        I=I2
        ELSE
        E(I2)=0.
        I=I1
        ENDIF
        P(1,I)=P(1,I1)+P(1,I2)
        P(2,I)=P(2,I1)+P(2,I2)
        P(3,I)=P(3,I1)+P(3,I2)
        scheck=(E10+E20)**2-p1**2-p2**2-p3**2
        if(scheck.lt.0) then
           write(99,*) 'scheck17: ', scheck
           write(99,*) 'scheck17', scheck,E10,E20,P(1,I),P(2,I),P(3,I)
           write(99,*) 'scheck17-1',E(I1),P(1,I1),P(2,I1),P(3,I1)
           write(99,*) 'scheck17-2',E(I2),P(1,I2),P(2,I2),P(3,I2)
           scheck=0.d0
        endif
        DM=SQRT(sngl(scheck))
        E(I)=DM
        RETURN
        END
