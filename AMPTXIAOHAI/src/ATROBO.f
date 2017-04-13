        SUBROUTINE ATROBO(THE,PHI,BEX,BEY,BEZ,IMIN,IMAX,IERROR)
        COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)
        DIMENSION ROT(3,3),PV(3)
        DOUBLE PRECISION DP(4),DBEX,DBEY,DBEZ,DGA,DGA2,DBEP,DGABEP
        SAVE   
        IERROR=0
              IF(IMIN.LE.0 .OR. IMAX.GT.N .OR. IMIN.GT.IMAX) RETURN
              IF(THE**2+PHI**2.GT.1E-20) THEN
           ROT(1,1)=COS(THE)*COS(PHI)
           ROT(1,2)=-SIN(PHI)
           ROT(1,3)=SIN(THE)*COS(PHI)
           ROT(2,1)=COS(THE)*SIN(PHI)
           ROT(2,2)=COS(PHI)
           ROT(2,3)=SIN(THE)*SIN(PHI)
           ROT(3,1)=-SIN(THE)
           ROT(3,2)=0.
           ROT(3,3)=COS(THE)
           DO 120 I=IMIN,IMAX
              DO 100 J=1,3
 100                 PV(J)=P(I,J)
                 DO 110 J=1,3
 110                    P(I,J)=ROT(J,1)*PV(1)+ROT(J,2)*PV(2)
     &                     +ROT(J,3)*PV(3)
 120                 CONTINUE
        ENDIF
              IF(BEX**2+BEY**2+BEZ**2.GT.1E-20) THEN
                DBEX=dble(BEX)
                DBEY=dble(BEY)
                DBEZ=dble(BEZ)
                DGA2=1D0-DBEX**2-DBEY**2-DBEZ**2
                IF(DGA2.LE.0D0) THEN
                        IERROR=1
                        RETURN
                ENDIF
                DGA=1D0/DSQRT(DGA2)
                DO 140 I=IMIN,IMAX
                   DO 130 J=1,4
 130                  DP(J)=dble(P(I,J))
                   DBEP=DBEX*DP(1)+DBEY*DP(2)+DBEZ*DP(3)
                   DGABEP=DGA*(DGA*DBEP/(1D0+DGA)+DP(4))
                   P(I,1)=sngl(DP(1)+DGABEP*DBEX)
                   P(I,2)=sngl(DP(2)+DGABEP*DBEY)
                   P(I,3)=sngl(DP(3)+DGABEP*DBEZ)
                   P(I,4)=sngl(DGA*(DP(4)+DBEP))
140                   CONTINUE
              ENDIF
              RETURN
              END
