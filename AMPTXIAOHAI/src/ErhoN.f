       real function ErhoN(em1,em2,lb1,lb2,srt)
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
      data arrayb/0.15,0.20,0.05,0.175,0.025,0.125,0.1,0.20,
     &0.53,0.34,0.05,0.07,0.15,0.45,0.45,0.058,
     &0.08,0.12,0.08/
       pi=3.1415926
       xs=0
       do 1001 ir=1,19
       IF(IR.LE.8)THEN
       if( ((lb1*lb2.eq.27.AND.(LB1.EQ.1.OR.LB2.EQ.1)).OR.
     &     (LB1*LB2.EQ.25*2.AND.(LB1.EQ.2.OR.LB2.EQ.2)))
     &       .OR.((lb1*lb2.eq.-25.AND.(LB1.EQ.-1.OR.LB2.EQ.-1)).OR.
     &     (LB1*LB2.EQ.-27*2.AND.(LB1.EQ.-2.OR.LB2.EQ.-2))) )
     &     branch=0.
        if((iabs(lb1*lb2).eq.26.AND.(iabs(LB1).EQ.1.OR.iabs(LB2).EQ.1))
     &   .OR.(iabs(LB1*LB2).EQ.26*2
     &   .AND.(iabs(LB1).EQ.2.OR.iabs(LB2).EQ.2)))
     &     branch=1./3.
       if( ((lb1*lb2.eq.27*2.AND.(LB1.EQ.2.OR.LB2.EQ.2)).OR.
     &     (LB1*LB2.EQ.25.AND.(LB1.EQ.1.OR.LB2.EQ.1)))
     &  .OR.((lb1*lb2.eq.-25*2.AND.(LB1.EQ.-2.OR.LB2.EQ.-2)).OR.
     &     (LB1*LB2.EQ.-27.AND.(LB1.EQ.-1.OR.LB2.EQ.-1))) )
     &     branch=2./3.
       ELSE
       if( ((lb1*lb2.eq.27.AND.(LB1.EQ.1.OR.LB2.EQ.1)).OR.
     &     (LB1*LB2.EQ.25*2.AND.(LB1.EQ.2.OR.LB2.EQ.2)))
     &       .OR.((lb1*lb2.eq.-25.AND.(LB1.EQ.-1.OR.LB2.EQ.-1)).OR.
     &     (LB1*LB2.EQ.-27*2.AND.(LB1.EQ.-2.OR.LB2.EQ.-2))) )
     &     branch=1.
        if((iabs(lb1*lb2).eq.26.AND.(iabs(LB1).EQ.1.OR.iabs(LB2).EQ.1))
     &   .OR.(iabs(LB1*LB2).EQ.26*2
     &   .AND.(iabs(LB1).EQ.2.OR.iabs(LB2).EQ.2)))
     &     branch=2./3.
       if( ((lb1*lb2.eq.27*2.AND.(LB1.EQ.2.OR.LB2.EQ.2)).OR.
     &     (LB1*LB2.EQ.25.AND.(LB1.EQ.1.OR.LB2.EQ.1)))
     &  .OR.((lb1*lb2.eq.-25*2.AND.(LB1.EQ.-2.OR.LB2.EQ.-2)).OR.
     &     (LB1*LB2.EQ.-27.AND.(LB1.EQ.-1.OR.LB2.EQ.-1))) )
     &     branch=1./3.
       ENDIF
       xs0=fdR(arraym(ir),arrayj(ir),arrayl(ir),
     &arrayw(ir),arrayb(ir),srt,EM1,EM2)
       xs=xs+1.3*pi*branch*xs0*(0.1973)**2
 1001 continue
       Erhon=xs
       return
       end
