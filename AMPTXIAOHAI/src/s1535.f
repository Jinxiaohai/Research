       real function s1535(SRT)
      SAVE   
       S0=2.424
       s1535=0.
       IF(SRT.LE.S0)RETURN
       S1535=2.*0.102*(SRT-S0)/(0.058+(SRT-S0)**2)
       return
       end
