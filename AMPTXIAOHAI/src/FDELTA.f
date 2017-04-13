        REAL FUNCTION FDELTA(DMASS)
      SAVE   
        AMN=0.938869
        AVPI=0.13803333
        AM0=1.232
        FD=0.25*WIDTH(DMASS)**2/((DMASS-1.232)**2
     1  +0.25*WIDTH(DMASS)**2)
        FDELTA=FD
        RETURN
        END
