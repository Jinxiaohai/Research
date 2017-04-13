      subroutine sdbelastic(SRT,sdb)
      PARAMETER (srt0=2.012)
      common/leadng/lb1,px1,py1,pz1,em1,e1,xfnl,yfnl,zfnl,tfnl,
     1     px1n,py1n,pz1n,dp1n
      common /dpi/em2,lb2
      common /para8/ idpert,npertd,idxsec
      SAVE   
      sdb=0.
      sdbel=0.
      if(srt.le.(em1+em2)) return
      s=srt**2
      if(idxsec.eq.1.or.idxsec.eq.3) then
         sdbel=fdbel(s)
      elseif(idxsec.eq.2.or.idxsec.eq.4) then
         threshold=em1+em2
         snew=(srt-threshold+srt0)**2
         sdbel=fdbel(snew)
      endif
      sdb=sdbel
      return
      end
