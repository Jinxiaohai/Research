        SUBROUTINE srree(lb1,lb2,srt)
        parameter (ETAM=0.5475,arho=0.77)
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
      SAVE   
        rree=0.
        if(lb1.ge.25.and.lb1.le.27.and.
     1       lb2.ge.25.and.lb2.le.27) then
           if(srt.gt.(2*ETAM)) rree=rrtoee(srt)
        elseif(lb1.eq.0.and.lb2.eq.0) then
           if(srt.gt.(2*arho)) rree=eetorr(srt)
        endif
        return
        END
