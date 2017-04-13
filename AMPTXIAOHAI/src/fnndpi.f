      function fnndpi(s)
      parameter(srt0=2.012)
      if(s.le.srt0**2) then
         fnndpi=0.
      else
         fnndpi=26.*exp(-(s-4.65)**2/0.1)+4.*exp(-(s-4.65)**2/2.)
     1        +0.28*exp(-(s-6.)**2/10.)
      endif
      return
      end
