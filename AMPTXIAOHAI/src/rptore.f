      real function rptore(srt)
      parameter (pimass=0.140,ETAM=0.5475,arho=0.77)
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
      SAVE   
      s2=srt**2
      pf2=(s2-(arho+ETAM)**2)*(s2-(arho-ETAM)**2)/2/sqrt(s2)
      pi2=(s2-(arho+pimass)**2)*(s2-(arho-pimass)**2)/2/sqrt(s2)
      rptore=1./3.*pf2/pi2*retorp(srt)
      return
      END
