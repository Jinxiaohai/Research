      real function xop2oe(srt)
      parameter (pimass=0.140,ETAM=0.5475,aomega=0.782)
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
      SAVE   
      s2=srt**2
      pf2=(s2-(aomega+ETAM)**2)*(s2-(aomega-ETAM)**2)/2/sqrt(s2)
      pi2=(s2-(aomega+pimass)**2)*(s2-(aomega-pimass)**2)/2/sqrt(s2)
      xop2oe=1./3.*pf2/pi2*xoe2op(srt)
      return
      END
