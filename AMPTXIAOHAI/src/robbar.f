      real function robbar(srt)
      parameter (pimass=0.140,arho=0.77,aomega=0.782)
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
      SAVE   
      sppb5p=xppbar(srt)*factr2(5)*ene**3/fsum
      pi2=(s-(arho+aomega)**2)*(s-(arho-aomega)**2)/4/s
      robbar=4./27.*sppb5p/pi2*wtot
      return
      END
