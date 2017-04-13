      real function pp2(srt)
      SAVE   
           pmass=0.9383 
       scheck=((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2
       if(scheck.lt.0) then
          write(99,*) 'scheck33: ', scheck
          scheck=0.
       endif
       plab=sqrt(scheck)
       pmin=2.
       pmax=2050
       if(plab.gt.pmax)then
       pp2=8.
       return
       endif
        if(plab .lt. pmin)then
        pp2 = 25.
        return
        end if
       a=11.2
       b=25.5
       c=0.151
       d=-1.62
       an=-1.12
        pp2 = a+b*(plab**an)+c*(alog(plab))**2+d*alog(plab)
       if(pp2.le.0)pp2=0
        return
        END
