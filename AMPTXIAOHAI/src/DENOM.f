        REAL FUNCTION DENOM(SRT,CON)
        PARAMETER (AP1=0.13496,
     1  AP2=0.13957,PI=3.1415926,AVMASS=0.9383)
      SAVE   
        AVPI=(AP1+2.*AP2)/3.
        AM0=1.232
        AMN=AVMASS
        AMP=AVPI
        AMAX=SRT-AVMASS
        AMIN=AVMASS+AVPI
        NMAX=200
        DMASS=(AMAX-AMIN)/FLOAT(NMAX)
        SUM=0.
        DO 10 I=1,NMAX+1
        DM=AMIN+FLOAT(I-1)*DMASS
        IF(CON.EQ.1.)THEN
        Q2=((DM**2-AMN**2+AMP**2)/(2.*DM))**2-AMP**2
           IF(Q2.GT.0.)THEN
           Q=SQRT(Q2)
           ELSE
           Q=1.E-06
           ENDIF
        TQ=0.47*(Q**3)/(AMP**2*(1.+0.6*(Q/AMP)**2))
        ELSE if(con.eq.2)then
        TQ=0.2
        AM0=1.44
       else if(con.eq.-1.)then
       tq=0.1
       am0=1.535
        ENDIF
        A1=4.*TQ*AM0**2/(AM0**2*TQ**2+(DM**2-AM0**2)**2)
        S=SRT**2
        P0=(S+DM**2-AMN**2)**2/(4.*S)-DM**2
        IF(P0.LE.0.)THEN
        P1=1.E-06
        ELSE
        P1=SQRT(P0)
        ENDIF
        F=DM*A1*P1
        IF((I.EQ.1).OR.(I.EQ.(NMAX+1)))THEN
        SUM=SUM+F*0.5
        ELSE
        SUM=SUM+F
        ENDIF
10      CONTINUE
        DENOM=SUM*DMASS/(2.*PI)
        RETURN
        END
