        REAL FUNCTION Fmassr(DMASS)
      SAVE   
        AM0=0.77
       wid=0.153
        Fmassr=am0*Wid/((DMASS**2-am0**2)**2
     1  +am0**2*Wid**2)
        RETURN
        END
