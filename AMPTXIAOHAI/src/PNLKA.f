      real function PNLKA(srt)
      SAVE   
      ala=1.116
      aka=0.498
      ana=0.939
      t1=ala+aka      
      if(srt.le.t1) THEN
      Pnlka=0
      Else
      IF(SRT.LT.1.7)sbbk=(0.9/0.091)*(SRT-T1)
      IF(SRT.GE.1.7)sbbk=0.09/(SRT-1.6)
      Pnlka=0.25*sbbk
       pnlka=pnlka/10.
      endif     
      return
      end
