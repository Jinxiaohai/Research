        SUBROUTINE ATTRAD(IERROR)
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
        COMMON/HIJDAT/HIDAT0(10,10),HIDAT(10)
        COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)
      COMMON/RNDF77/NSEED
        SAVE   
        IERROR=0
40        SM=0.
        JL=1
        DO 30 I=1,N-1
           S=2.*(P(I,4)*P(I+1,4)-P(I,1)*P(I+1,1)-P(I,2)*P(I+1,2)
     &                -P(I,3)*P(I+1,3))+P(I,5)**2+P(I+1,5)**2
           IF(S.LT.0.) S=0.
           WP=SQRT(S)-1.5*(P(I,5)+P(I+1,5))
           IF(WP.GT.SM) THEN
              PBT1=P(I,1)+P(I+1,1)
              PBT2=P(I,2)+P(I+1,2)
              PBT3=P(I,3)+P(I+1,3)
              PBT4=P(I,4)+P(I+1,4)
              BTT=(PBT1**2+PBT2**2+PBT3**2)/PBT4**2
              IF(BTT.GE.1.0-1.0E-10) GO TO 30
              IF((I.NE.1.OR.I.NE.N-1).AND.
     &             (K(I,2).NE.21.AND.K(I+1,2).NE.21)) GO TO 30
              JL=I
              SM=WP
           ENDIF
30        CONTINUE
        S=(SM+1.5*(P(JL,5)+P(JL+1,5)))**2
              IF(SM.LT.HIPR1(5)) GOTO 2
              IF(JL+1.EQ.N) GOTO 190
              DO 160 J=N,JL+2,-1
                      K(J+1,1)=K(J,1)
                K(J+1,2)=K(J,2)
                      DO 150 M=1,5
150                           P(J+1,M)=P(J,M)
160                   CONTINUE
190           N=N+1
              P1=P(JL,1)+P(JL+1,1)
              P2=P(JL,2)+P(JL+1,2)
              P3=P(JL,3)+P(JL+1,3)
              P4=P(JL,4)+P(JL+1,4)
              BEX=-P1/P4
              BEY=-P2/P4
              BEZ=-P3/P4
        IMIN=JL
        IMAX=JL+1
              CALL ATROBO(0.,0.,BEX,BEY,BEZ,IMIN,IMAX,IERROR)
        IF(IERROR.NE.0) RETURN
              CTH=P(JL,3)/SQRT(P(JL,4)**2-P(JL,5)**2)
              IF(ABS(CTH).GT.1.0)  CTH=MAX(-1.,MIN(1.,CTH))
              THETA=ACOS(CTH)
              PHI=ULANGL(P(JL,1),P(JL,2))
              CALL ATROBO(0.,-PHI,0.,0.,0.,IMIN,IMAX,IERROR)
              CALL ATROBO(-THETA,0.,0.,0.,0.,IMIN,IMAX,IERROR)
1        CALL AR3JET(S,X1,X3,JL)
              CALL ARORIE(S,X1,X3,JL)                
        IF(HIDAT(2).GT.0.0) THEN
                 PTG1=SQRT(P(JL,1)**2+P(JL,2)**2)
                 PTG2=SQRT(P(JL+1,1)**2+P(JL+1,2)**2)
                 PTG3=SQRT(P(JL+2,1)**2+P(JL+2,2)**2)
           PTG=MAX(PTG1,PTG2,PTG3)
           IF(PTG.GT.HIDAT(2)) THEN
              FMFACT=EXP(-(PTG**2-HIDAT(2)**2)/HIPR1(2)**2)
              IF(RANART(NSEED).GT.FMFACT) GO TO 1
           ENDIF
        ENDIF
        IMIN=JL
        IMAX=JL+2
              CALL ATROBO(THETA,PHI,-BEX,-BEY,-BEZ,IMIN,IMAX,IERROR)
        IF(IERROR.NE.0) RETURN
              K(JL+2,1)=K(JL+1,1)
        K(JL+2,2)=K(JL+1,2)
        K(JL+2,3)=K(JL+1,3)
        K(JL+2,4)=K(JL+1,4)
        K(JL+2,5)=K(JL+1,5)
              P(JL+2,5)=P(JL+1,5)
              K(JL+1,1)=2
        K(JL+1,2)=21
        K(JL+1,3)=0
        K(JL+1,4)=0
        K(JL+1,5)=0
              P(JL+1,5)=0.
             IF(SM.GE.HIPR1(5)) GOTO 40
2              K(1,1)=2
        K(1,3)=0
        K(1,4)=0
        K(1,5)=0
              K(N,1)=1
        K(N,3)=0
        K(N,4)=0
        K(N,5)=0
              RETURN
              END
