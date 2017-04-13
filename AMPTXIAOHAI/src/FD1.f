      REAL FUNCTION FD1(DMASS,aj,al,width,widb0,srt)
      SAVE   
        AMN=0.938
        AmP=0.138
       amd=amn
           Ak02= 0.25*(DMASS**2-amd**2-amp**2)**2
     &           -(Amp*amd)**2
            IF (ak02 .GT. 0.) THEN
              Q0 = SQRT(ak02/DMASS)
            ELSE
              Q0= 0.0
             fd1=0
           return
            END IF
           Ak2= 0.25*(srt**2-amd**2-amp**2)**2
     &           -(Amp*amd)**2
            IF (ak2 .GT. 0.) THEN
              Q = SQRT(ak2/DMASS)
            ELSE
              Q= 0.00
             fd1=0
             return
            END IF
       b=widb0*1.2*dmass/srt*(q/q0)**(2.*al+1)
     &  /(1.+0.2*(q/q0)**(2*al))
        FD1=(2.*aj+1)*WIDTH**2*b/((srt-dmass)**2
     1  +0.25*WIDTH**2)/(2.*q**2)
        RETURN
        END
