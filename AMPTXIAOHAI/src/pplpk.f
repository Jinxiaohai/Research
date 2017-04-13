      real function pplpk(srt)
      SAVE   
           pmass=0.9383 
       pplpk=0.
       scheck=((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2
       if(scheck.lt.0) then
          write(99,*) 'scheck35: ', scheck
          scheck=0.
       endif
       plab=sqrt(scheck)
       pmin=2.82
       pmax=25.0
       if(plab.gt.pmax)then
       pplpk=0.036
       return
       endif
        if(plab .lt. pmin)then
        pplpk = 0.
        return
        end if
       a=0.0654
       b=-3.16
       c=-0.0029
       an=-4.14
        pplpk = a+b*(plab**an)+c*(alog(plab))**2
       if(pplpk.le.0)pplpk=0
        return
        END
