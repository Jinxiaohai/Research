       subroutine ppxS(lb1,lb2,srt,ppsig,spprho,ipp)
       parameter (amp=0.14,pi=3.1415926)
      SAVE   
       PPSIG=0.0
        spprho=0.0
       IPP=0
       IF(SRT.LE.0.3)RETURN
       q=sqrt((srt/2)**2-amp**2)
       esigma=5.8*amp
       tsigma=2.06*q
       erho=0.77
       trho=0.095*q*(q/amp/(1.+(q/erho)**2))**2
       esi=esigma-srt
       if(esi.eq.0)then
       d00=pi/2.
       go to 10
       endif
       d00=atan(tsigma/2./esi)
10       erh=erho-srt
       if(erh.eq.0.)then
       d11=pi/2.
       go to 20
       endif
       d11=atan(trho/2./erh)
20       d20=-0.12*q/amp
       s0=8.*pi*sin(d00)**2/q**2
       s1=8*pi*3*sin(d11)**2/q**2
       s2=8*pi*5*sin(d20)**2/q**2
        s0=s0*0.197**2*10.
        s1=s1*0.197**2*10.
        s2=s2*0.197**2*10.
       spprho=s1/2.
       IF(LB1.EQ.5.AND.LB2.EQ.5)THEN
       IPP=1
       PPSIG=S2
       RETURN
       ENDIF
       IF((LB1.EQ.5.AND.LB2.EQ.4).OR.(LB1.EQ.4.AND.LB2.EQ.5))THEN
       IPP=2
       PPSIG=S2/2.+S1/2.
       RETURN
       ENDIF
       IF((LB1.EQ.5.AND.LB2.EQ.3).OR.(LB1.EQ.3.AND.LB2.EQ.5))THEN
       IPP=3
       PPSIG=S2/6.+S1/2.+S0/3.
       RETURN
       ENDIF
       IF(LB1.EQ.4.AND.LB2.EQ.4)THEN
       IPP=4
       PPSIG=2*S2/3.+S0/3.
       RETURN
       ENDIF
       IF((LB1.EQ.4.AND.LB2.EQ.3).OR.(LB1.EQ.3.AND.LB2.EQ.4))THEN
       IPP=5
       PPSIG=S2/2.+S1/2.
       RETURN
       ENDIF
       IF(LB1.EQ.3.AND.LB2.EQ.3)THEN
       IPP=6
       PPSIG=S2
       ENDIF
       return
       end
