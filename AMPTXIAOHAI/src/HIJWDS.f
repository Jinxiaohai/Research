        SUBROUTINE HIJWDS(IA,IDH,XHIGH)
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
        COMMON/WOOD/R,D,FNORM,W
        DIMENSION IAA(20),RR(20),DD(20),WW(20)
        EXTERNAL RWDSAX,WDSAX
        SAVE   
        DATA IAA/2,4,12,16,27,32,40,56,63,93,184,197,208,7*0./
        DATA RR/0.01,.964,2.355,2.608,2.84,3.458,3.766,3.971,4.214,
     1        4.87,6.51,6.38,6.624,7*0./
        DATA DD/0.5882,.322,.522,.513,.569,.61,.586,.5935,.586,.573,
     1        .535,.535,.549,7*0./
        DATA WW/0.0,.517,-0.149,-0.051,0.,-0.208,-0.161,13*0./
              A=IA
              D=0.54
        R=1.19*A**(1./3.) - 1.61*A**(-1./3.)
        W=0.
        DO 10 I=1,13
                IF (IA.EQ.IAA(I)) THEN
                        R=RR(I)
                             D=DD(I)
                              W=WW(I)
                      END IF
10            CONTINUE
              FNORM=1.0
              XLOW=0.
              XHIGH=R+ 12.*D
              IF (W.LT.-0.01)  THEN
                      IF (XHIGH.GT.R/SQRT(ABS(W))) XHIGH=R/SQRT(ABS(W))
              END IF
              FGAUS=GAUSS1(RWDSAX,XLOW,XHIGH,0.001)
              FNORM=1./FGAUS
        IF (IDH.EQ.1) THEN
           HINT1(72)=R
           HINT1(73)=D
           HINT1(74)=W
           HINT1(75)=FNORM/4.0/HIPR1(40)
        ELSE IF (IDH.EQ.2) THEN
           HINT1(76)=R
           HINT1(77)=D
           HINT1(78)=W
           HINT1(79)=FNORM/4.0/HIPR1(40)
        ENDIF
              CALL HIFUN(IDH,XLOW,XHIGH,RWDSAX)
              RETURN
              END
