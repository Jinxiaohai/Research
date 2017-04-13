        REAL FUNCTION FDE(DMASS,SRT,CON)
      SAVE   
        AMN=0.938869
        AVPI=0.13803333
        AM0=1.232
        FD=4.*(AM0**2)*WIDTH(DMASS)/((DMASS**2-1.232**2)**2
     1  +AM0**2*WIDTH(DMASS)**2)
        IF(CON.EQ.1.)THEN
        P11=(SRT**2+DMASS**2-AMN**2)**2
     1  /(4.*SRT**2)-DMASS**2
       if(p11.le.0)p11=1.E-06
       p1=sqrt(p11)
        ELSE
        DMASS=AMN+AVPI
        P11=(SRT**2+DMASS**2-AMN**2)**2
     1  /(4.*SRT**2)-DMASS**2
       if(p11.le.0)p11=1.E-06
       p1=sqrt(p11)
        ENDIF
        FDE=FD*P1*DMASS
        RETURN
        END
