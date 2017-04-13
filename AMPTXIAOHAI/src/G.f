        FUNCTION G(Y1,Y2,PT2)
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
        CALL PARTON(F,X1,X2,PT2)
        G11=( (F(1,1)+F(1,2))*(F(2,3)+F(2,4)+F(2,5)+F(2,6))
     &      +(F(1,3)+F(1,4))*(F(2,5)+F(2,6)) )*SUBCR1(T,U)
        G12=( (F(2,1)+F(2,2))*(F(1,3)+F(1,4)+F(1,5)+F(1,6))
     &      +(F(2,3)+F(2,4))*(F(1,5)+F(1,6)) )*SUBCR1(U,T)
        G13=(F(1,1)*F(2,1)+F(1,2)*F(2,2)+F(1,3)*F(2,3)+F(1,4)*F(2,4)
     &      +F(1,5)*F(2,5)+F(1,6)*F(2,6))*(SUBCR1(U,T)
     &      +SUBCR1(T,U)-8.D0/T/U/27.D0)
        G2=(AF-1)*(F(1,1)*F(2,2)+F(2,1)*F(1,2)+F(1,3)*F(2,4)
     &     +F(2,3)*F(1,4)+F(1,5)*F(2,6)+F(2,5)*F(1,6))*SUBCR2(T,U)
        G31=(F(1,1)*F(2,2)+F(1,3)*F(2,4)+F(1,5)*F(2,6))*SUBCR3(T,U)
        G32=(F(2,1)*F(1,2)+F(2,3)*F(1,4)+F(2,5)*F(1,6))*SUBCR3(U,T)
        G4=(F(1,1)*F(2,2)+F(2,1)*F(1,2)+F(1,3)*F(2,4)+F(2,3)*F(1,4)+
     1        F(1,5)*F(2,6)+F(2,5)*F(1,6))*SUBCR4(T,U)
        G5=AF*F(1,7)*F(2,7)*SUBCR5(T,U)
        G61=F(1,7)*(F(2,1)+F(2,2)+F(2,3)+F(2,4)+F(2,5)
     &      +F(2,6))*SUBCR6(T,U)
        G62=F(2,7)*(F(1,1)+F(1,2)+F(1,3)+F(1,4)+F(1,5)
     &      +F(1,6))*SUBCR6(U,T)
        G7=F(1,7)*F(2,7)*SUBCR7(T,U)
        G=(G11+G12+G13+G2+G31+G32+G4+G5+G61+G62+G7)*dble(HIPR1(17))*
     1        3.14159D0*APH**2/SS**2
        RETURN
        END
