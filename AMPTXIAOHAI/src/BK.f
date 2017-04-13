        FUNCTION BK(X)
        COMMON /BESEL/X4
        SAVE   
        BK=EXP(-X)*(X**2-X4**2)**2.50/15.0
        RETURN
        END
