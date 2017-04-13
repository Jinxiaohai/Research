        REAL FUNCTION SIGMA(SRT,ID,IOI,IOF)
        PARAMETER (AMU=0.9383,AMP=0.1384,PI=3.1415926,HC=0.19733)
      SAVE   
        IF(ID.EQ.1)THEN
        AMASS0=1.22
        T0 =0.12
        ELSE
        AMASS0=1.43
        T0 =0.2
        ENDIF
        IF((IOI.EQ.1).AND.(IOF.EQ.1))THEN
        ALFA=3.772
        BETA=1.262
        AM0=1.188
        T=0.09902
        ENDIF
        IF((IOI.EQ.1).AND.(IOF.EQ.0))THEN
        ALFA=15.28
        BETA=0.
        AM0=1.245
        T=0.1374
        ENDIF
        IF((IOI.EQ.0).AND.(IOF.EQ.1))THEN
        ALFA=146.3
        BETA=0.
        AM0=1.472
        T=0.02649
        ENDIF
        ZPLUS=(SRT-AMU-AMASS0)*2./T0
        ZMINUS=(AMU+AMP-AMASS0)*2./T0
        deln=ATAN(ZPLUS)-ATAN(ZMINUS)
       if(deln.eq.0)deln=1.E-06
        AMASS=AMASS0+(T0/4.)*ALOG((1.+ZPLUS**2)/(1.+ZMINUS**2))
     1  /deln
        S=SRT**2
        P2=S/4.-AMU**2
        S0=(AMU+AM0)**2
        P02=S0/4.-AMU**2
        P0=SQRT(P02)
        PR2=(S-(AMU-AMASS)**2)*(S-(AMU+AMASS)**2)/(4.*S)
        IF(PR2.GT.1.E-06)THEN
        PR=SQRT(PR2)
        ELSE
        PR=0.
        SIGMA=1.E-06
        RETURN
        ENDIF
        SS=AMASS**2
        Q2=(SS-(AMU-AMP)**2)*(SS-(AMU+AMP)**2)/(4.*SS)
        IF(Q2.GT.1.E-06)THEN
        Q=SQRT(Q2)
        ELSE
        Q=0.
        SIGMA=1.E-06
        RETURN
        ENDIF
        SS0=AM0**2
        Q02=(SS0-(AMU-AMP)**2)*(SS0-(AMU+AMP)**2)/(4.*SS0)
        scheck=Q02
        if(scheck.lt.0) then
           write(99,*) 'scheck20: ', scheck
           scheck=0.
        endif
        Q0=SQRT(scheck)
        SIGMA=PI*(HC)**2/(2.*P2)*ALFA*(PR/P0)**BETA*AM0**2*T**2
     1  *(Q/Q0)**3/((SS-AM0**2)**2+AM0**2*T**2)
        SIGMA=SIGMA*10.
       IF(SIGMA.EQ.0)SIGMA=1.E-06
        RETURN
        END
