      real function ppk0(srt)
      real   xarray(7), earray(7)
      SAVE   
      data xarray /0.030,0.025,0.025,0.026,0.02,0.014,0.06/
      data earray /3.67,4.95,5.52,6.05,6.92,7.87,10./
           pmass=0.9383 
       ppk0=0
       if(srt.le.2.63)return
       if(srt.gt.4.54)then
       ppk0=0.037
       return
       endif
        plab=sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)
        if (plab .lt. earray(1)) then
        ppk0 = xarray(1)
        return
      end if
      do 1001 ie = 1,7
        if (earray(ie) .eq. plab) then
          ppk0 = xarray(ie)
          go to 10
        else if (earray(ie) .gt. plab) then
          ymin = alog(xarray(ie-1))
          ymax = alog(xarray(ie))
          xmin = alog(earray(ie-1))
          xmax = alog(earray(ie))
          ppk0 = exp(ymin + (alog(plab)-xmin)*(ymax-ymin)
     &/(xmax-xmin) )
          go to 10
        end if
 1001 continue
10       continue
      return
        END
