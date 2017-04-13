        real function akNPsg(pkaon)
      SAVE   
         if(pkaon.le.0.345)then
           sigma1=0.624*pkaon**(-1.83)
         else
           sigma1=0.7*pkaon**(-2.09)
         endif
        akNPsg=sigma1
        return
        end   
