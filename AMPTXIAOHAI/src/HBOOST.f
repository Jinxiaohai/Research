        SUBROUTINE HBOOST
              IMPLICIT DOUBLE PRECISION(D)  
              COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5) 
              COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
        SAVE   
        DO 100 I=1,N
           DBETA=dble(P(I,3)/P(I,4))
           IF(ABS(DBETA).GE.1.D0) THEN
              DB=dble(HINT1(2))
              IF(DB.GT.0.99999999D0) THEN 
                 WRITE(6,*) '(HIBOOT:) boost vector too large' 
                 DB=0.99999999D0
              ENDIF 
              DGA=1D0/SQRT(1D0-DB**2)
              DP3=dble(P(I,3))
              DP4=dble(P(I,4))
              P(I,3)=sngl((DP3+DB*DP4)*DGA)
              P(I,4)=sngl((DP4+DB*DP3)*DGA)
              GO TO 100
           ENDIF
           Y=0.5*sngl(DLOG((1.D0+DBETA)/(1.D0-DBETA)))
           AMT=SQRT(P(I,1)**2+P(I,2)**2+P(I,5)**2)
           P(I,3)=AMT*SINH(Y+HINT1(3))
           P(I,4)=AMT*COSH(Y+HINT1(3))
100        CONTINUE
        RETURN
        END
