        FUNCTION SUBCR7(T,U)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        SUBCR7=9.D0/2.D0*(3.D0-T*U-U/T**2-T/U**2)
        RETURN
        END
