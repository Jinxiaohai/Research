      SUBROUTINE crdbel(PX,PY,PZ,SRT,I1,I2,IBLOCK,
     1     NTAG,sig,NT,ianti)
      PARAMETER (MAXSTR=150001,MAXR=1)
      COMMON /AA/R(3,MAXSTR)
      COMMON /BB/ P(3,MAXSTR)
      COMMON /BG/BETAX,BETAY,BETAZ,GAMMA
      COMMON /CC/ E(MAXSTR)
      COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
      COMMON /AREVT/ IAEVT, IARUN, MISS
      common/leadng/lb1,px1,py1,pz1,em1,e1,xfnl,yfnl,zfnl,tfnl,
     1     px1n,py1n,pz1n,dp1n
      common /dpi/em2,lb2
      common /para8/ idpert,npertd,idxsec
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
      SAVE   
      IBLOCK=0
      NTAG=0
      EM1=E(I1)
      EM2=E(I2)
      s=srt**2
      if(sig.le.0) return
      IBLOCK=503
      if(iabs(lb1).eq.42) then
         ideut=i1
         lbb=lb2
         idb=i2
      else
         ideut=i2
         lbb=lb1
         idb=i1
      endif
      if((idpert.eq.1.or.idpert.eq.2).and.dpertp(ideut).ne.1.) then
         if(ianti.eq.0) then
            write(91,*) '  d+',lbb,' (pert d B elastic) @nt=',nt
     1           ,' @prob=',dpertp(ideut),p(1,idb),p(2,idb)
     2           ,p(1,ideut),p(2,ideut)
         else
            write(91,*) '  d+',lbb,' (pert dbar Bbar elastic) @nt=',nt
     1           ,' @prob=',dpertp(ideut),p(1,idb),p(2,idb)
     2           ,p(1,ideut),p(2,ideut)
         endif
         scheck=(s-(em1+em2)**2)*(s-(em1-em2)**2)
         if(scheck.lt.0) then
            write(99,*) 'scheck53: ', scheck
            scheck=0.
         endif
         pfinal=sqrt(scheck)/2./srt
         CALL dbelangle(pxn,pyn,pzn,pfinal)
         CALL ROTATE(PX,PY,PZ,Pxn,Pyn,Pzn)
         EdCM=SQRT(E(ideut)**2+Pxn**2+Pyn**2+Pzn**2)
         PdBETA=Pxn*BETAX+Pyn*BETAY+Pzn*BETAZ
         TRANSF=GAMMA*(GAMMA*PdBETA/(GAMMA+1.)+EdCM)
         Pt1d=BETAX*TRANSF+Pxn
         Pt2d=BETAY*TRANSF+Pyn
         Pt3d=BETAZ*TRANSF+Pzn
         p(1,ideut)=pt1d
         p(2,ideut)=pt2d
         p(3,ideut)=pt3d
         PX1=P(1,I1)
         PY1=P(2,I1)
         PZ1=P(3,I1)
         ID(I1)=2
         ID(I2)=2
         R(1,ideut)=R(1,idb)
         R(2,ideut)=R(2,idb)
         R(3,ideut)=R(3,idb)
         return
      endif
      if(ianti.eq.0) then
         write (91,*) ' d+',lbb,' (regular d B elastic) @evt#',
     1        iaevt,' @nt=',nt,' lb1,2=',lb1,lb2
      else
         write (91,*) ' d+',lbb,' (regular dbar Bbar elastic) @evt#',
     1        iaevt,' @nt=',nt,' lb1,2=',lb1,lb2
      endif
      scheck=(s-(em1+em2)**2)*(s-(em1-em2)**2)
      if(scheck.lt.0) then
         write(99,*) 'scheck54: ', scheck
         scheck=0.
      endif
      pfinal=sqrt(scheck)/2./srt
      CALL dbelangle(pxn,pyn,pzn,pfinal)
      CALL ROTATE(PX,PY,PZ,Pxn,Pyn,Pzn)
      E1CM=SQRT(E(I1)**2+Pxn**2+Pyn**2+Pzn**2)
      P1BETA=Pxn*BETAX+Pyn*BETAY+Pzn*BETAZ
      TRANSF=GAMMA*(GAMMA*P1BETA/(GAMMA+1.)+E1CM)
      Pt1i1=BETAX*TRANSF+Pxn
      Pt2i1=BETAY*TRANSF+Pyn
      Pt3i1=BETAZ*TRANSF+Pzn
      p(1,i1)=pt1i1
      p(2,i1)=pt2i1
      p(3,i1)=pt3i1
      E2CM=SQRT(E(I2)**2+Pxn**2+Pyn**2+Pzn**2)
      P2BETA=-Pxn*BETAX-Pyn*BETAY-Pzn*BETAZ
      TRANSF=GAMMA*(GAMMA*P2BETA/(GAMMA+1.)+E2CM)
      Pt1I2=BETAX*TRANSF-Pxn
      Pt2I2=BETAY*TRANSF-Pyn
      Pt3I2=BETAZ*TRANSF-Pzn
      p(1,i2)=pt1i2
      p(2,i2)=pt2i2
      p(3,i2)=pt3i2
      PX1=P(1,I1)
      PY1=P(2,I1)
      PZ1=P(3,I1)
      EM1=E(I1)
      EM2=E(I2)
      ID(I1)=2
      ID(I2)=2
      RETURN
      END
