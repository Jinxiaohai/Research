        FUNCTION SUBCR4(T,U)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        SUBCR4=8.D0/3.D0*(T**2+U**2)*(4.D0/9.D0/T/U-1.D0)
        RETURN
        END
