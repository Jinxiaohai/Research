      real function xppbar(srt)
       Parameter (pmass=0.9383,xmax=400.)
      SAVE   
       xppbar=1.e-06
       plab2=(srt**2/(2.*pmass)-pmass)**2-pmass**2
       if(plab2.gt.0)then
           plab=sqrt(plab2)
       xppbar=67./(plab**0.7)
       if(xppbar.gt.xmax)xppbar=xmax
       endif
         return
      END
