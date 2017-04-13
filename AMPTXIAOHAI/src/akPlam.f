        real function akPlam(pkaon)
      SAVE   
        p=pkaon
        if(pkaon.lt.0.2.or. pkaon.ge.10.0) sigma=0.
        if(pkaon.ge.0.2.and.pkaon.lt.0.9) sigma=50.*p**2-67.*p+24.
        if(pkaon.ge.0.9.and.pkaon.lt.10.0) sigma=3.0*pkaon**(-2.6)
        akPlam=sigma
        return
        end
