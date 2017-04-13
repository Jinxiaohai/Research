        REAL FUNCTION FD5(DMASS,SRT,CON)
      SAVE   
        AMN=0.938869
        AVPI=0.13803333
        AM0=1.535
        FD=4.*(AM0**2)*W1535(DMASS)/((DMASS**2-1.535**2)**2
     1  +AM0**2*W1535(DMASS)**2)
        IF(CON.EQ.1.)THEN
           scheck=(SRT**2+DMASS**2-AMN**2)**2/(4.*SRT**2)-DMASS**2
           if(scheck.lt.0) then
              write(99,*) 'scheck11: ', scheck
              scheck=0.
           endif
           P1=SQRT(scheck)
        ELSE
        DMASS=AMN+AVPI
        scheck=(SRT**2+DMASS**2-AMN**2)**2/(4.*SRT**2)-DMASS**2
        if(scheck.lt.0) then
           write(99,*) 'scheck12: ', scheck
           scheck=0.
        endif
        P1=SQRT(scheck)
        ENDIF
        FD5=FD*P1*DMASS
        RETURN
        END
