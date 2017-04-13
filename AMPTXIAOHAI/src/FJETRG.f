        FUNCTION FJETRG(X,WGT)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        REAL HIPR1(100),HINT1(100),PTMAX,PTMIN
        COMMON/HPARNT/HIPR1,IHPR2(50),HINT1,IHNT2(50)
        DIMENSION X(10)
        SAVE   
        PTMIN=ABS(HIPR1(10))-0.25
        PTMIN=MAX(PTMIN,HIPR1(8))
        AM2=0.D0
        IF(IHPR2(3).EQ.3) THEN
           AM2=dble(HIPR1(7)**2)
           PTMIN=MAX(0.0,HIPR1(10))
        ENDIF
        PTMAX=ABS(HIPR1(10))+0.25
        IF(HIPR1(10).LE.0.0) PTMAX=HINT1(1)/2.0-sngl(AM2)
        IF(PTMAX.LE.PTMIN) PTMAX=PTMIN+0.25
        PT2=dble(PTMAX**2-PTMIN**2)*X(1)+dble(PTMIN)**2
        AMT2=PT2+AM2
        XT=2.0d0*DSQRT(AMT2)/dble(HINT1(1))
        YMX1=DLOG(1.0d0/XT+DSQRT(1.0d0/XT**2-1.0d0))
        Y1=2.0d0*YMX1*X(2)-YMX1
        YMX2=DLOG(2.0d0/XT-DEXP(Y1))
        YMN2=DLOG(2.0d0/XT-DEXP(-Y1))
        Y2=(YMX2+YMN2)*X(3)-YMN2
        IF(IHPR2(3).EQ.3) THEN
           GTRIG=2.0d0*GHVQ(Y1,Y2,AMT2)
        ELSE IF(IHPR2(3).EQ.2) THEN
           GTRIG=2.0d0*GPHOTN(Y1,Y2,PT2)
        ELSE
           GTRIG=G(Y1,Y2,PT2)
        ENDIF
        FJETRG=2.0d0*YMX1*(YMX2+YMN2)*dble(PTMAX**2-PTMIN**2)
     &                *GTRIG/2.0d0
        RETURN
        END
