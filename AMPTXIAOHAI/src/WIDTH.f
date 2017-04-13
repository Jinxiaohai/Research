        REAL FUNCTION WIDTH(DMASS)
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
            WIDTH = 0.47 * QAVAIL**3 /
     &              (PIMASS**2 * (1.+0.6*(QAVAIL/PIMASS)**2))
        RETURN
        END
