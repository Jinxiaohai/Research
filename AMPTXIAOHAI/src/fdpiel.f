      function fdpiel(s)
      parameter(srt0=2.012)
      if(s.le.srt0**2) then
         fdpiel=0.
      else
         fdpiel=63.*exp(-(s-4.67)**2/0.15)+15.*exp(-(s-6.25)**2/0.3)
      endif
      return
      end
