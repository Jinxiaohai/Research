        REAL FUNCTION FRHO(DMASS)
      SAVE   
        AM0=0.77
       WID=0.153
        FD=0.25*wid**2/((DMASS-AM0)**2+0.25*WID**2)
        FRHO=FD
        RETURN
        END
