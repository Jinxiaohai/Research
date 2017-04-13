       Subroutine M1535(LB1,LB2,SRT,X1535)
      SAVE   
       S0=2.424
       x1535=0.
       IF(SRT.LE.S0)RETURN
       SIGMA=2.*0.102*(SRT-S0)/(0.058+(SRT-S0)**2)
       IF((LB1*LB2.EQ.18.AND.(LB1.EQ.2.OR.LB2.EQ.2)).OR.
     &     (LB1*LB2.EQ.6.AND.(LB1.EQ.1.OR.LB2.EQ.1)).or.
     &     (lb1*lb2.eq.8.AND.(LB1.EQ.1.OR.LB2.EQ.1)))then
       X1535=SIGMA
       return
       ENDIF
       IF(LB1*LB2.EQ.7)THEN
       X1535=3.*SIGMA
       RETURN
       ENDIF 
       IF((LB1*LB2.EQ.11).OR.
     &     (LB1*LB2.EQ.20.AND.(LB1.EQ.2.OR.LB2.EQ.2)))THEN
       X1535=SIGMA
       RETURN
       ENDIF
       IF((LB1*LB2.EQ.10.AND.(LB1.EQ.1.OR.LB2.EQ.1)).OR.
     &     (LB1*LB2.EQ.22.AND.(LB1.EQ.2.OR.LB2.EQ.2)))
     &     X1535=3.*SIGMA
       RETURN
       END
