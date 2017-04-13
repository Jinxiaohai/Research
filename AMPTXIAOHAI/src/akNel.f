        real function akNel(pkaon)
      SAVE   
        if(pkaon.lt.0.5.or. pkaon.ge.4.0) sigma1=0.
        if(pkaon.ge.0.5.and.pkaon.lt.1.0) sigma1=20.*pkaon**2.74
        if(pkaon.ge.1.0.and.pkaon.lt.4.0) sigma1=20.*pkaon**(-1.8)
        akNel=sigma1
        return
        end
