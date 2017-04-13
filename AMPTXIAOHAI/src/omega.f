      real function omega(srt)
      SAVE   
       pmass=0.9383
       omass=0.782
       tomega=0.0084
       omega=0.00000001
       if(srt.le.2.68)return
       ESMIN=2.*0.9383+omass-tomega/2.
       es=srt
       omega=0.36*(es-esmin)/(1.25+(es-esmin)**2)
       return
       end
