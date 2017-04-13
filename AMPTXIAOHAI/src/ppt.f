      real function ppt(srt)
      SAVE   
           pmass=0.9383 
       scheck=((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2
       if(scheck.lt.0) then
          write(99,*) 'scheck34: ', scheck
          scheck=0.
       endif
       plab=sqrt(scheck)
       pmin=3. 
       pmax=2100
      if ((plab .lt. pmin).or.(plab.gt.pmax)) then
        ppt = 55.
        return
      end if
       a=45.6
       b=219.0
       c=0.410
       d=-3.41
       an=-4.23
        ppt = a+b*(plab**an)+c*(alog(plab))**2+d*alog(plab)
       if(ppt.le.0)ppt=0.0
        return
        END
