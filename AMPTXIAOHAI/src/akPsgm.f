        real function akPsgm(pkaon)
      SAVE   
        if(pkaon.lt.0.2.or. pkaon.ge.1.5) sigma1=0.
        if(pkaon.ge.0.2.and.pkaon.lt.1.5) sigma1=0.6*pkaon**(-1.8)
        akPsgm=sigma1
        return
        end
