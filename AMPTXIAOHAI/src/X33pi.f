      real function X33pi(srt)
      real   xarray(12), earray(12)
      SAVE   
      data xarray /0.02,0.22,0.74,1.10,1.76,1.84,2.20,
     &2.40,2.15,2.60,2.30,1.70/
      data earray /2.23,2.81,3.67,4.00,4.95,5.52,5.97,
     &6.04,6.60,6.90,10.01,19./
           pmass=0.9383 
       x33pi=1.E-06
       if(srt.le.2.3)return
      plab=sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)
      if (plab .lt. earray(1)) then
        x33pi = xarray(1)
        return
      end if
      do 1001 ie = 1,12
        if (earray(ie) .eq. plab) then
          x33pi= xarray(ie)
          return
        else if (earray(ie) .gt. plab) then
          ymin = alog(xarray(ie-1))
          ymax = alog(xarray(ie))
          xmin = alog(earray(ie-1))
          xmax = alog(earray(ie))
          x33pi= exp(ymin + (alog(plab)-xmin)*(ymax-ymin)
     &    /(xmax-xmin))
          return
        end if
 1001   continue
        return
        END
