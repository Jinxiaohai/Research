        FUNCTION GHVQ(Y1,Y2,AMT2)
        IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
        REAL HIPR1(100),HINT1(100)
        COMMON/HPARNT/HIPR1,IHPR2(50),HINT1,IHNT2(50)
        DIMENSION F(2,7)
        SAVE   
        XT=2.0d0*DSQRT(AMT2)/dble(HINT1(1))
        X1=0.5d0*XT*(DEXP(Y1)+DEXP(Y2))
        X2=0.5d0*XT*(DEXP(-Y1)+DEXP(-Y2))
        SS=X1*X2*dble(HINT1(1))**2
        AF=4.0d0
        IF(IHPR2(18).NE.0) AF=5.0d0
        DLAM=dble(HIPR1(15))
        APH=12.0d0*3.1415926d0/(33.d0-2.d0*AF)/DLOG(AMT2/DLAM**2)
        CALL PARTON(F,X1,X2,AMT2)
        Gqq=4.d0*(DCOSH(Y1-Y2)+dble(HIPR1(7))**2/AMT2)
     &       /(1.D0+DCOSH(Y1-Y2))
     &       /9.d0*(F(1,1)*F(2,2)+F(1,2)*F(2,1)+F(1,3)*F(2,4)
     &       +F(1,4)*F(2,3)+F(1,5)*F(2,6)+F(1,6)*F(2,5))
        Ggg=(8.D0*DCOSH(Y1-Y2)-1.D0)
     &       *(DCOSH(Y1-Y2)+2.d0*dble(HIPR1(7))**2
     &       /AMT2-2.d0*dble(HIPR1(7))**4/AMT2**2)/(1.d0+DCOSH(Y1-Y2))
     &       /24.d0*F(1,7)*F(2,7)
        GHVQ=(Gqq+Ggg)*dble(HIPR1(23))*3.14159d0*APH**2/SS**2
        RETURN
        END
