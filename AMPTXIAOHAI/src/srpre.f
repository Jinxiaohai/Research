        SUBROUTINE srpre(lb1,lb2,srt)
        parameter (pimass=0.140,ETAM=0.5475,arho=0.77)
        common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
        common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
      SAVE   
        rpre=0.
        if(lb1.ge.25.and.lb1.le.27.and.lb2.ge.3.and.lb2.le.5) then
           if(srt.gt.(ETAM+arho)) rpre=rptore(srt)
        elseif(lb2.ge.25.and.lb2.le.27.and.lb1.ge.3.and.lb1.le.5) then
           if(srt.gt.(ETAM+arho)) rpre=rptore(srt)
        elseif(lb1.ge.25.and.lb1.le.27.and.lb2.eq.0) then
           if(srt.gt.(pimass+arho)) rpre=retorp(srt)
        elseif(lb2.ge.25.and.lb2.le.27.and.lb1.eq.0) then
           if(srt.gt.(pimass+arho)) rpre=retorp(srt)
        endif
        return
        END
