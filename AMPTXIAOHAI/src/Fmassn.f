        REAL FUNCTION Fmassn(DMASS)
      SAVE   
        AM0=1.44
        Fmassn=am0*W1440(DMASS)/((DMASS**2-am0**2)**2
     1  +am0**2*W1440(DMASS)**2)
        RETURN
        END
