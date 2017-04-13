      real function X3pi(srt)
      real   xarray(12), earray(12)
      SAVE   
      data xarray /0.02,0.4,1.15,1.60,2.19,2.85,2.30,
     &3.10,2.47,2.60,2.40,1.70/
      data earray /2.23,2.81,3.67,4.00,4.95,5.52,5.97,
     &6.04,6.60,6.90,10.01,19./
           pmass=0.9383 
       x3pi=1.E-06
       if(srt.le.2.3)return
      plab=sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)
      if (plab .lt. earray(1)) then
        x3pi = xarray(1)
        return
      end if
      do 1001 ie = 1,12
        if (earray(ie) .eq. plab) then
          x3pi= xarray(ie)
          return
        else if (earray(ie) .gt. plab) then
          ymin = alog(xarray(ie-1))
          ymax = alog(xarray(ie))
          xmin = alog(earray(ie-1))
          xmax = alog(earray(ie))
          X3pi= exp(ymin + (alog(plab)-xmin)*(ymax-ymin)
     &                                            /(xmax-xmin) )
          return
        end if
 1001 continue
      return
        END
