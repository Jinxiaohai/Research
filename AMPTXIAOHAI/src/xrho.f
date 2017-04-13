      real function xrho(srt)
      SAVE   
       pmass=0.9383
       rmass=0.77
       trho=0.151
       xrho=0.000000001
       if(srt.le.2.67)return
       ESMIN=2.*0.9383+rmass-trho/2.
       ES=srt
       xrho0=0.24*(es-esmin)/(1.4+(es-esmin)**2)
       xrho=3.*Xrho0
       return
       end
