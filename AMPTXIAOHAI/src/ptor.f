      real function ptor(srt)
      parameter (pimass=0.140,arho=0.77)
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
      SAVE   
      s2=srt**2
      ptor=9*(s2-4*arho**2)/(s2-4*pimass**2)*rtop(srt)
      return
      END
