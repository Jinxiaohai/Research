        REAL FUNCTION W1535(DMASS)
      SAVE   
        AVMASS=0.938868
        PIMASS=0.137265
           AUX = 0.25*(DMASS**2-AVMASS**2-PIMASS**2)**2
     &           -(AVMASS*PIMASS)**2
            IF (AUX .GT. 0.) THEN
              QAVAIL = SQRT(AUX / DMASS**2)
            ELSE
              QAVAIL = 1.E-06
            END IF
            W1535 = 0.15* QAVAIL/0.467
        RETURN
        END
