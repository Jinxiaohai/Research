        REAL FUNCTION FNS(DMASS,SRT,CON)
      SAVE   
        WIDTH=0.2
        AMN=0.938869
        AVPI=0.13803333
        AN0=1.43
        FN=4.*(AN0**2)*WIDTH/((DMASS**2-1.44**2)**2+AN0**2*WIDTH**2)
        IF(CON.EQ.1.)THEN
           scheck=(SRT**2+DMASS**2-AMN**2)**2/(4.*SRT**2)-DMASS**2
           if(scheck.lt.0) then
              write(99,*) 'scheck13: ', scheck
              scheck=0.
           endif
           P1=SQRT(scheck)
        ELSE
        DMASS=AMN+AVPI
        scheck=(SRT**2+DMASS**2-AMN**2)**2/(4.*SRT**2)-DMASS**2
        if(scheck.lt.0) then
           write(99,*) 'scheck14: ', scheck
           scheck=0.
        endif
        P1=SQRT(scheck)
        ENDIF
        FNS=FN*P1*DMASS
        RETURN
        END
