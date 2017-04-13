        FUNCTION SUBCR3(T,U)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        SUBCR3=4.D0/9.D0*(T**2+U**2+(1.D0+U**2)/T**2
     1        -2.D0*U**2/3.D0/T)
        RETURN
        END
