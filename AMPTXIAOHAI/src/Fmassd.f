        REAL FUNCTION Fmassd(DMASS)
      SAVE   
        AM0=1.232
        Fmassd=am0*WIDTH(DMASS)/((DMASS**2-am0**2)**2
     1  +am0**2*WIDTH(DMASS)**2)
        RETURN
        END
