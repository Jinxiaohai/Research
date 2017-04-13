        SUBROUTINE spppe(lb1,lb2,srt)
        parameter (pimass=0.140,ETAM=0.5475)
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
      SAVE   
        pppe=0.
        if((lb1.ge.3.and.lb1.le.5).and.(lb2.ge.3.and.lb2.le.5)) then
           if(srt.gt.(ETAM+pimass)) pppe=pptope(srt)
        elseif((lb1.ge.3.and.lb1.le.5).and.lb2.eq.0) then
           pppe=petopp(srt)
        elseif((lb2.ge.3.and.lb2.le.5).and.lb1.eq.0) then
           pppe=petopp(srt)
        endif
        return
        END
