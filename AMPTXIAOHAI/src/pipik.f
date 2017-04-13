      real function pipik(srt)
      real   xarray(5), earray(5)
      SAVE   
      data xarray /0.001, 0.7,1.5,1.7,2.0/
      data earray /1.,1.2,1.6,2.0,2.4/
           pmass=0.9383 
       pipik=0.
       if(srt.le.1.)return
       if(srt.gt.2.4)then
           pipik=2.0/2.
           return
       endif
        if (srt .lt. earray(1)) then
           pipik =xarray(1)/2.
           return
        end if
      do 1001 ie = 1,5
        if (earray(ie) .eq. srt) then
          pipik = xarray(ie)
          go to 10
        else if (earray(ie) .gt. srt) then
          ymin = alog(xarray(ie-1))
          ymax = alog(xarray(ie))
          xmin = alog(earray(ie-1))
          xmax = alog(earray(ie))
          pipik = exp(ymin + (alog(srt)-xmin)*(ymax-ymin)
     &/(xmax-xmin) )
          go to 10
        end if
 1001 continue
10       PIPIK=PIPIK/2.
       continue
      return
        END
