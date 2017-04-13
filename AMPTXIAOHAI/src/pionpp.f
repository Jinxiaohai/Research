      real function pionpp(srt)
      SAVE   
           pmass=0.14 
       pmass1=0.938
       PIONPP=0.00001
       IF(SRT.LE.1.22)RETURN
        plab=sqrt(((srt**2-pmass**2-pmass1**2)/(2.*pmass1))**2-pmass**2)
       pmin=0.3
       pmax=25.0
       if(plab.gt.pmax)then
       pionpp=20./10.
       return
       endif
        if(plab .lt. pmin)then
        pionpp = 0.
        return
        end if
       a=24.3
       b=-12.3
       c=0.324
       an=-1.91
       d=-2.44
        pionpp = a+b*(plab**an)+c*(alog(plab))**2+d*alog(plab)
       if(pionpp.le.0)pionpp=0
       pionpp=pionpp/10.
        return
        END
