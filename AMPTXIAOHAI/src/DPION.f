      REAL FUNCTION DPION(EM1,EM2,LB1,LB2,SRT)
      dimension   arrayj(19),arrayl(19),arraym(19),
     &arrayw(19),arrayb(19)
      SAVE   
      data arrayj /0.5,1.5,0.5,0.5,2.5,2.5,1.5,0.5,1.5,3.5,
     &1.5,0.5,1.5,0.5,2.5,0.5,1.5,2.5,3.5/
      data arrayl/1,2,0,0,2,3,2,1,1,3,
     &1,0,2,0,3,1,1,2,3/
      data arraym /1.44,1.52,1.535,1.65,1.675,1.68,1.70,1.71,
     &1.72,1.99,1.60,1.62,1.70,1.90,1.905,1.910,
     &1.86,1.93,1.95/
      data arrayw/0.2,0.125,0.15,0.15,0.155,0.125,0.1,0.11,
     &0.2,0.29,0.25,0.16,0.28,0.15,0.3,0.22,0.25,
     &0.25,0.24/
      data arrayb/0.15,0.25,0.,0.05,0.575,0.125,0.379,0.10,
     &0.10,0.062,0.45,0.60,0.6984,0.05,0.25,0.089,
     &0.19,0.2,0.13/
       pi=3.1415926
       amn=0.94
       amp=0.14
       xs=0
       do 1001 ir=1,19
       BRANCH=0.
       if(ir.LE.8)THEN
       IF( ((LB1*LB2.EQ.5*7.AND.(LB1.EQ.5.OR.LB2.EQ.5)).OR.
     &     (LB1*LB2.EQ.3*8.AND.(LB1.EQ.3.OR.LB2.EQ.3)))
     &       .OR.((LB1*LB2.EQ.-3*7.AND.(LB1.EQ.3.OR.LB2.EQ.3)).OR.
     &     (LB1*LB2.EQ.-5*8.AND.(LB1.EQ.5.OR.LB2.EQ.5))) )
     &     branch=1./6.
       IF((iabs(LB1*LB2).EQ.4*7.AND.(LB1.EQ.4.OR.LB2.EQ.4)).OR.
     &     (iabs(LB1*LB2).EQ.4*8.AND.(LB1.EQ.4.OR.LB2.EQ.4)))
     &     branch=1./3.
       IF( ((LB1*LB2.EQ.5*6.AND.(LB1.EQ.5.OR.LB2.EQ.5)).OR.
     &     (LB1*LB2.EQ.3*9.AND.(LB1.EQ.3.OR.LB2.EQ.3)))
     &       .OR.((LB1*LB2.EQ.-3*6.AND.(LB1.EQ.3.OR.LB2.EQ.3)).OR.
     &     (LB1*LB2.EQ.-5*9.AND.(LB1.EQ.5.OR.LB2.EQ.5))) )
     &     branch=1./2.
       ELSE
       IF( ((LB1*LB2.EQ.5*8.AND.(LB1.EQ.5.OR.LB2.EQ.5)).OR.
     &     (LB1*LB2.EQ.5*6.AND.(LB1.EQ.5.OR.LB2.EQ.5)))
     &        .OR.((LB1*LB2.EQ.-3*8.AND.(LB1.EQ.3.OR.LB2.EQ.3)).OR.
     &     (LB1*LB2.EQ.-3*6.AND.(LB1.EQ.3.OR.LB2.EQ.3))) )
     &     branch=2./5.
       IF( ((LB1*LB2.EQ.3*9.AND.(LB1.EQ.3.OR.LB2.EQ.3)).OR.
     &     (LB1*LB2.EQ.3*7.AND.(LB1.EQ.3.OR.LB2.EQ.3)))
     &        .OR. ((LB1*LB2.EQ.-5*9.AND.(LB1.EQ.5.OR.LB2.EQ.5)).OR.
     &     (LB1*LB2.EQ.-5*7.AND.(LB1.EQ.5.OR.LB2.EQ.5))) )
     &     branch=2./5.
       IF( ((LB1*LB2.EQ.5*7.AND.(LB1.EQ.5.OR.LB2.EQ.5)).OR.
     &     (LB1*LB2.EQ.3*8.AND.(LB1.EQ.3.OR.LB2.EQ.3)))
     &        .OR.((LB1*LB2.EQ.-3*7.AND.(LB1.EQ.3.OR.LB2.EQ.3)).OR.
     &     (LB1*LB2.EQ.-5*8.AND.(LB1.EQ.5.OR.LB2.EQ.5))) )
     &     branch=8./15.
       IF((iabs(LB1*LB2).EQ.4*7.AND.(LB1.EQ.4.OR.LB2.EQ.4)).OR.
     &     (iabs(LB1*LB2).EQ.4*8.AND.(LB1.EQ.4.OR.LB2.EQ.4)))
     &     branch=1./15.
       IF((iabs(LB1*LB2).EQ.4*9.AND.(LB1.EQ.4.OR.LB2.EQ.4)).OR.
     &     (iabs(LB1*LB2).EQ.4*6.AND.(LB1.EQ.4.OR.LB2.EQ.4)))
     &     branch=3./5.
       ENDIF
       xs0=fd2(arraym(ir),arrayj(ir),arrayl(ir),
     &arrayw(ir),arrayb(ir),EM1,EM2,srt)
       xs=xs+1.3*pi*branch*xs0*(0.1973)**2
 1001   continue
       DPION=XS
       RETURN
       end
