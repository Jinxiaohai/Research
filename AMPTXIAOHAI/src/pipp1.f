      real function pipp1(srt)
      SAVE   
           pmass=0.14 
       pmass1=0.938
       PIPP1=0.0001
       IF(SRT.LE.1.22)RETURN
        plab=sqrt(((srt**2-pmass**2-pmass1**2)/(2.*pmass1))**2-pmass**2)
       pmin=0.3
       pmax=25.0
       if(plab.gt.pmax)then
       pipp1=20./10.
       return
       endif
        if(plab .lt. pmin)then
        pipp1 = 0.
        return
        end if
       a=26.6
       b=-7.18
       c=0.327
       an=-1.86
       d=-2.81
        pipp1 = a+b*(plab**an)+c*(alog(plab))**2+d*alog(plab)
       if(pipp1.le.0)pipp1=0
       PIPP1=PIPP1/10.
        return
        END
