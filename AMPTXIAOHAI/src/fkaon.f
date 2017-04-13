       Real function fkaon(p,pmax)
      SAVE   
       fmax=0.148
       if(pmax.eq.0.)pmax=0.000001
       fkaon=(1.-p/pmax)*(p/pmax)**2
       if(fkaon.gt.fmax)fkaon=fmax
       fkaon=fkaon/fmax
       return
       end
