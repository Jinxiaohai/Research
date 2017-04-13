      function fdbel(s)
      parameter(srt0=2.012)
      if(s.le.srt0**2) then
         fdbel=0.
      else
         fdbel=2500.*exp(-(s-7.93)**2/0.003)
     1        +300.*exp(-(s-7.93)**2/0.1)+10.
      endif
      return
      end
