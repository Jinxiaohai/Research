      real function eetorr(srt)
      parameter (ETAM=0.5475,arho=0.77)
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
      SAVE   
      s2=srt**2
      eetorr=81.*(s2-4*arho**2)/(s2-4*etam**2)*rrtoee(srt)
      return
      END
