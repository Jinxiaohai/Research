      real function FOURPI(srt)
      real   xarray(10), earray(10)
      SAVE   
      data xarray /0.0001,1.986597,6.411932,7.636956,    
     &9.598362,9.889740,10.24317,10.80138,11.86988,12.83925/    
      data earray /2.468,2.718,2.968,0.322E+01,
     &0.347E+01, 0.372E+01, 0.397E+01, 0.422E+01, 0.447E+01,
     &0.472E+01/
           pmass=0.14 
       pmass1=0.938
       FOURPI=0.000001
       if(srt.le.1.52)return
        plab=SRT
      if (plab .lt. earray(1)) then
        FOURPI= 0.00001
        return
      end if
      do 1001 ie = 1,10
        if (earray(ie) .eq. plab) then
          FOURPI= xarray(ie)
          return
        else if (earray(ie) .gt. plab) then
          ymin = alog(xarray(ie-1))
          ymax = alog(xarray(ie))
          xmin = alog(earray(ie-1))
          xmax = alog(earray(ie))
          FOURPI= exp(ymin + (alog(plab)-xmin)*(ymax-ymin)
     &    /(xmax-xmin) )
          return
        end if
 1001   continue
      return
        END
