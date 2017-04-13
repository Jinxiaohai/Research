      real function TWOPI(srt)
      real   xarray(19), earray(19)
      SAVE   
      data xarray /0.300E-05,0.187E+01,0.110E+02,0.149E+02,0.935E+01,
     &0.765E+01,0.462E+01,0.345E+01,0.241E+01,0.185E+01,0.165E+01,
     &0.150E+01,0.132E+01,0.117E+01,0.116E+01,0.100E+01,0.856E+00,
     &0.745E+00,0.300E-05/
      data earray /0.122E+01, 0.147E+01, 0.172E+01, 0.197E+01,
     &0.222E+01, 0.247E+01, 0.272E+01, 0.297E+01, 0.322E+01,
     &0.347E+01, 0.372E+01, 0.397E+01, 0.422E+01, 0.447E+01,
     &0.472E+01, 0.497E+01, 0.522E+01, 0.547E+01, 0.572E+01/
           pmass=0.14 
       pmass1=0.938
       TWOPI=0.000001
       if(srt.le.1.22)return
        plab=SRT
      if (plab .lt. earray(1)) then
        TWOPI= 0.00001
        return
      end if
      do 1001 ie = 1,19
        if (earray(ie) .eq. plab) then
          TWOPI= xarray(ie)
          return
        else if (earray(ie) .gt. plab) then
          ymin = alog(xarray(ie-1))
          ymax = alog(xarray(ie))
          xmin = alog(earray(ie-1))
          xmax = alog(earray(ie))
          TWOPI= exp(ymin + (alog(plab)-xmin)*(ymax-ymin)
     &    /(xmax-xmin) )
          return
        end if
 1001   continue
      return
        END
