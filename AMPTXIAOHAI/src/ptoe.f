      real function ptoe(srt)
      parameter (pimass=0.140,ETAM=0.5475)
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
      SAVE   
      s2=srt**2
      ptoe=1./9.*(s2-4*etam**2)/(s2-4*pimass**2)*etop(srt)
      return
      END
