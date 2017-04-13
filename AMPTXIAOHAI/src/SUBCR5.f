        FUNCTION SUBCR5(T,U)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        SUBCR5=3.D0/8.D0*(T**2+U**2)*(4.D0/9.D0/T/U-1.D0)
        RETURN
        END
