        SUBROUTINE sopoe(lb1,lb2,srt)
        parameter (ETAM=0.5475,aomega=0.782)
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
      SAVE   
        xopoe=0.
        if((lb1.eq.28.and.lb2.ge.3.and.lb2.le.5).or.
     1       (lb2.eq.28.and.lb1.ge.3.and.lb1.le.5)) then
           if(srt.gt.(aomega+ETAM)) xopoe=xop2oe(srt)
        elseif((lb1.eq.28.and.lb2.eq.0).or.
     1          (lb1.eq.0.and.lb2.eq.28)) then
           if(srt.gt.(aomega+ETAM)) xopoe=xoe2op(srt)
        endif
        return
        END
