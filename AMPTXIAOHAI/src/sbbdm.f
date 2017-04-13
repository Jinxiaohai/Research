      subroutine sbbdm(srt,sdprod,ianti,lbm,xmm,pfinal)
      PARAMETER (xmd=1.8756,AP1=0.13496,AP2=0.13957,
     1     xmrho=0.770,xmomega=0.782,xmeta=0.548,srt0=2.012)
      common/leadng/lb1,px1,py1,pz1,em1,e1,xfnl,yfnl,zfnl,tfnl,
     1     px1n,py1n,pz1n,dp1n
      common /dpi/em2,lb2
      common /para8/ idpert,npertd,idxsec
      COMMON/RNDF77/NSEED
      SAVE   
      sdprod=0.
      sbbdpi=0.
      sbbdrho=0.
      sbbdomega=0.
      sbbdeta=0.
      if(srt.le.(em1+em2)) return
      ilb1=iabs(lb1)
      ilb2=iabs(lb2)
      s=srt**2
      scheck=(s-(em1+em2)**2)*(s-(em1-em2)**2)
      if(scheck.le.0) then
         write(99,*) 'scheck50: ', scheck
         stop
      endif
      pinitial=sqrt(scheck)/2./srt
      fs=fnndpi(s)
      if(idxsec.eq.1.or.idxsec.eq.2) then
      else
         if(ilb1.ge.1.and.ilb1.le.2.and.
     1        ilb2.ge.1.and.ilb2.le.2) then
            pifactor=9./8.
         elseif((ilb1.ge.1.and.ilb1.le.2.and.
     1           ilb2.ge.6.and.ilb2.le.9).or.
     2           (ilb2.ge.1.and.ilb2.le.2.and.
     1           ilb1.ge.6.and.ilb1.le.9)) then
            pifactor=9./64.
         elseif((ilb1.ge.1.and.ilb1.le.2.and.
     1           ilb2.ge.10.and.ilb2.le.13).or.
     2           (ilb2.ge.1.and.ilb2.le.2.and.
     1           ilb1.ge.10.and.ilb1.le.13)) then
            pifactor=9./16.
         elseif(ilb1.ge.6.and.ilb1.le.9.and.
     1           ilb2.ge.6.and.ilb2.le.9) then
            pifactor=9./128.
         elseif((ilb1.ge.6.and.ilb1.le.9.and.
     1           ilb2.ge.10.and.ilb2.le.13).or.
     2           (ilb2.ge.6.and.ilb2.le.9.and.
     1           ilb1.ge.10.and.ilb1.le.13)) then
            pifactor=9./64.
         elseif((ilb1.ge.10.and.ilb1.le.11.and.
     1           ilb2.ge.10.and.ilb2.le.11).or.
     2           (ilb2.ge.12.and.ilb2.le.13.and.
     1           ilb1.ge.12.and.ilb1.le.13)) then
            pifactor=9./8.
         elseif((ilb1.ge.10.and.ilb1.le.11.and.
     1           ilb2.ge.12.and.ilb2.le.13).or.
     2           (ilb2.ge.10.and.ilb2.le.11.and.
     1           ilb1.ge.12.and.ilb1.le.13)) then
            pifactor=9./16.
         endif
      endif
      IF((ilb1*ilb2).EQ.1)THEN
         lbm=5
         if(ianti.eq.1) lbm=3
         xmm=ap2
      ELSEIF(ilb1.EQ.2.AND.ilb2.EQ.2)THEN
         lbm=3
         if(ianti.eq.1) lbm=5
         xmm=ap2
      ELSEIF((ilb1*ilb2).EQ.2)THEN
         lbm=4
         xmm=ap1
      ELSE
         lbm=3+int(3 * RANART(NSEED))
         if(lbm.eq.4) then
            xmm=ap1
         else
            xmm=ap2
         endif
      ENDIF
      if(srt.ge.(xmd+xmm)) then
         pfinal=sqrt((s-(xmd+xmm)**2)*(s-(xmd-xmm)**2))/2./srt
         if((ilb1.eq.1.and.ilb2.eq.1).or.
     1        (ilb1.eq.2.and.ilb2.eq.2)) then
            sbbdpi=fs*pfinal/pinitial/4.
         elseif((ilb1.eq.1.and.ilb2.eq.2).or.
     1           (ilb1.eq.2.and.ilb2.eq.1)) then
            sbbdpi=fs*pfinal/pinitial/4./2.
         else
            if(idxsec.eq.1) then
               sbbdpi=fs*pfinal/pinitial*3./16.
            elseif(idxsec.eq.2.or.idxsec.eq.4) then
               threshold=amax1(xmd+xmm,em1+em2)
               snew=(srt-threshold+srt0)**2
               if(idxsec.eq.2) then
                  sbbdpi=fnndpi(snew)*pfinal/pinitial*3./16.
               elseif(idxsec.eq.4) then
                  sbbdpi=fnndpi(snew)*pfinal/pinitial/6.*pifactor
               endif
            elseif(idxsec.eq.3) then
               sbbdpi=fs*pfinal/pinitial/6.*pifactor
            endif
         endif
      endif
      if(srt.gt.(xmd+xmrho)) then
         pfinal=sqrt((s-(xmd+xmrho)**2)*(s-(xmd-xmrho)**2))/2./srt
         if(idxsec.eq.1) then
            sbbdrho=fs*pfinal/pinitial*3./16.
         elseif(idxsec.eq.2.or.idxsec.eq.4) then
            threshold=amax1(xmd+xmrho,em1+em2)
            snew=(srt-threshold+srt0)**2
            if(idxsec.eq.2) then
               sbbdrho=fnndpi(snew)*pfinal/pinitial*3./16.
            elseif(idxsec.eq.4) then
               sbbdrho=fnndpi(snew)*pfinal/pinitial/6.*(pifactor*3.)
            endif
         elseif(idxsec.eq.3) then
            sbbdrho=fs*pfinal/pinitial/6.*(pifactor*3.)
         endif
      endif
      if(srt.gt.(xmd+xmomega)) then
         pfinal=sqrt((s-(xmd+xmomega)**2)*(s-(xmd-xmomega)**2))/2./srt
         if(idxsec.eq.1) then
            sbbdomega=fs*pfinal/pinitial*3./16.
         elseif(idxsec.eq.2.or.idxsec.eq.4) then
            threshold=amax1(xmd+xmomega,em1+em2)
            snew=(srt-threshold+srt0)**2
            if(idxsec.eq.2) then
               sbbdomega=fnndpi(snew)*pfinal/pinitial*3./16.
            elseif(idxsec.eq.4) then
               sbbdomega=fnndpi(snew)*pfinal/pinitial/6.*pifactor
            endif
         elseif(idxsec.eq.3) then
            sbbdomega=fs*pfinal/pinitial/6.*pifactor
         endif
      endif
      if(srt.gt.(xmd+xmeta)) then
         pfinal=sqrt((s-(xmd+xmeta)**2)*(s-(xmd-xmeta)**2))/2./srt
         if(idxsec.eq.1) then
            sbbdeta=fs*pfinal/pinitial*3./16.
         elseif(idxsec.eq.2.or.idxsec.eq.4) then
            threshold=amax1(xmd+xmeta,em1+em2)
            snew=(srt-threshold+srt0)**2
            if(idxsec.eq.2) then
               sbbdeta=fnndpi(snew)*pfinal/pinitial*3./16.
            elseif(idxsec.eq.4) then
               sbbdeta=fnndpi(snew)*pfinal/pinitial/6.*(pifactor/3.)
            endif
         elseif(idxsec.eq.3) then
            sbbdeta=fs*pfinal/pinitial/6.*(pifactor/3.)
         endif
      endif
      sdprod=sbbdpi+sbbdrho+sbbdomega+sbbdeta
      if(sdprod.le.0) return
      x1=RANART(NSEED)
      if(x1.le.sbbdpi/sdprod) then
      elseif(x1.le.(sbbdpi+sbbdrho)/sdprod) then
         lbm=25+int(3*RANART(NSEED))
         xmm=xmrho
      elseif(x1.le.(sbbdpi+sbbdrho+sbbdomega)/sdprod) then
         lbm=28
         xmm=xmomega
      else
         lbm=0
         xmm=xmeta
      endif
      return
      end
