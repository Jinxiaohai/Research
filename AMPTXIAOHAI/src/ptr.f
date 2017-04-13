       function ptr(ptmax,iseed)
        COMMON/TABLE/ xarray(0:1000),earray(0:1000)
      COMMON/RNDF77/NSEED
      SAVE   
       ptr=0.
       if(ptmax.le.1.e-02)then
       ptr=ptmax
       return
       endif
       if(ptmax.gt.2.01)ptmax=2.01
       tryial=ptdis(ptmax)/ptdis(2.01)
       XT=RANART(NSEED)*tryial
        do 50 ie = 1,200
        if (earray(ie) .eq. xT) then
          ptr = xarray(ie)
       return
       end if
          if(xarray(ie-1).le.0.00001)go to 50
          if(xarray(ie).le.0.00001)go to 50
          if(earray(ie-1).le.0.00001)go to 50
          if(earray(ie).le.0.00001)go to 50
        if (earray(ie) .gt. xT) then
          ymin = alog(xarray(ie-1))
          ymax = alog(xarray(ie))
          xmin = alog(earray(ie-1))
          xmax = alog(earray(ie))
          ptr= exp(ymin + (alog(xT)-xmin)*(ymax-ymin)
     &    /(xmax-xmin) )
       if(ptr.gt.ptmax)ptr=ptmax
       return
       endif
50      continue
       return
       end
