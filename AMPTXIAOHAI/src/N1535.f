       Subroutine N1535(LB1,LB2,SRT,X1535)
      SAVE   
       S0=2.424
       x1535=0.
       IF(SRT.LE.S0)RETURN
       SIGMA=2.*0.102*(SRT-S0)/(0.058+(SRT-S0)**2)
       IF((LB1*LB2.EQ.1).OR.
     &     (LB1.EQ.2.AND.LB2.EQ.2))then
       X1535=SIGMA
       return
       endif
       IF(LB1*LB2.EQ.2)then
       X1535=3.*SIGMA
       return
       endif 
       IF((LB1*LB2.EQ.63.AND.(LB1.EQ.7.OR.LB2.EQ.7)).OR.
     &     (LB1*LB2.EQ.64.AND.(LB1.EQ.8.OR.LB2.EQ.8)).OR.
     &     (LB1*LB2.EQ.48.AND.(LB1.EQ.6.OR.LB2.EQ.6)).OR.
     &     (LB1*LB2.EQ.49.AND.(LB1.EQ.7.OR.LB2.EQ.7)))then
       X1535=SIGMA
       return
       endif
       IF((LB1*LB2.EQ.54.AND.(LB1.EQ.6.OR.LB2.EQ.6)).OR.
     &     (LB1*LB2.EQ.56.AND.(LB1.EQ.7.OR.LB2.EQ.7)))then
       X1535=3.*SIGMA
       return
       endif
       IF((LB1.EQ.10.AND.LB2.EQ.10).OR.
     &     (LB1.EQ.11.AND.LB2.EQ.11))X1535=SIGMA
       IF(LB1*LB2.EQ.110.AND.(LB1.EQ.10.OR.LB2.EQ.10))X1535=3.*SIGMA
       RETURN
       END
