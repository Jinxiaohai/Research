        FUNCTION GPHOTN(Y1,Y2,PT2)
        IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
        REAL HIPR1(100),HINT1(100)
        COMMON/HPARNT/HIPR1,IHPR2(50),HINT1,IHNT2(50)
        DIMENSION F(2,7)
        SAVE   
        XT=2.d0*DSQRT(PT2)/dble(HINT1(1))
        X1=0.5d0*XT*(DEXP(Y1)+DEXP(Y2))
        X2=0.5d0*XT*(DEXP(-Y1)+DEXP(-Y2))
        Z=DSQRT(1.D0-XT**2/X1/X2)
        SS=X1*X2*dble(HINT1(1))**2
        T=-(1.d0-Z)/2.d0
        U=-(1.d0+Z)/2.d0
        AF=3.d0
        DLAM=dble(HIPR1(15))
        APH=12.d0*3.1415926d0/(33.d0-2.d0*AF)/DLOG(PT2/DLAM**2)
        APHEM=1.d0/137.d0
        CALL PARTON(F,X1,X2,PT2)
        G11=-(U**2+1.d0)/U/3.d0*F(1,7)*(4.d0*F(2,1)+4.d0*F(2,2)
     &      +F(2,3)+F(2,4)+F(2,5)+F(2,6))/9.d0
        G12=-(T**2+1.d0)/T/3.d0*F(2,7)*(4.d0*F(1,1)+4.d0*F(1,2)
     &      +F(1,3)+F(1,4)+F(1,5)+F(1,6))/9.d0
        G2=8.d0*(U**2+T**2)/U/T/9.d0*(4.d0*F(1,1)*F(2,2)
     &     +4.d0*F(1,2)*F(2,1)+F(1,3)*F(2,4)+F(1,4)*F(2,3)
     &     +F(1,5)*F(2,6)+F(1,6)*F(2,5))/9.d0
        GPHOTN=(G11+G12+G2)*dble(HIPR1(23))*3.14159d0*APH*APHEM/SS**2
        RETURN
        END
