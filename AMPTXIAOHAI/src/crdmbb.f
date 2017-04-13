      SUBROUTINE crdmbb(PX,PY,PZ,SRT,I1,I2,IBLOCK,
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
      common /dpifsl/lbnn1,lbnn2,lbnd1,lbnd2,lbns1,lbns2,lbnp1,lbnp2,
     1     lbdd1,lbdd2,lbds1,lbds2,lbdp1,lbdp2,lbss1,lbss2,
     2     lbsp1,lbsp2,lbpp1,lbpp2
      common /dpifsm/xmnn1,xmnn2,xmnd1,xmnd2,xmns1,xmns2,xmnp1,xmnp2,
     1     xmdd1,xmdd2,xmds1,xmds2,xmdp1,xmdp2,xmss1,xmss2,
     2     xmsp1,xmsp2,xmpp1,xmpp2
      common /dpisig/sdmel,sdmnn,sdmnd,sdmns,sdmnp,sdmdd,sdmds,sdmdp,
     1     sdmss,sdmsp,sdmpp
      COMMON/RNDF77/NSEED
      SAVE   
      IBLOCK=0
      NTAG=0
      EM1=E(I1)
      EM2=E(I2)
      s=srt**2
      if(sig.le.0) return
      if(iabs(lb1).eq.42) then
         ideut=i1
         lbm=lb2
         idm=i2
      else
         ideut=i2
         lbm=lb1
         idm=i1
      endif
      if((idpert.eq.1.or.idpert.eq.2).and.dpertp(ideut).ne.1.) then
         x1=RANART(NSEED)
         if(x1.le.sdmel/sig)then
            if(ianti.eq.0) then
               write(91,*) '  d+',lbm,' (pert d M elastic) @nt=',nt
     1              ,' @prob=',dpertp(ideut)
            else
               write(91,*) '  d+',lbm,' (pert dbar M elastic) @nt=',nt
     1              ,' @prob=',dpertp(ideut)
            endif
            scheck=(s-(em1+em2)**2)*(s-(em1-em2)**2)
            if(scheck.lt.0) then
               write(99,*) 'scheck51: ', scheck
               scheck=0.
            endif
            pfinal=sqrt(scheck)/2./srt
            CALL dmelangle(pxn,pyn,pzn,pfinal)
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
            IBLOCK=504
            PX1=P(1,I1)
            PY1=P(2,I1)
            PZ1=P(3,I1)
            ID(I1)=2
            ID(I2)=2
            R(1,ideut)=R(1,idm)
            R(2,ideut)=R(2,idm)
            R(3,ideut)=R(3,idm)
         else
            if(ianti.eq.0) then
               write(91,*) '  d+',lbm,' ->BB (pert d destrn) @nt=',nt
     1              ,' @prob=',dpertp(ideut)
            else
               write(91,*) '  d+',lbm,' ->BB (pert dbar destrn) @nt=',nt
     1              ,' @prob=',dpertp(ideut)
            endif
            e(ideut)=0.
            IBLOCK=502
         endif
         return
      endif
      IBLOCK=502
      x1=RANART(NSEED)
      if(x1.le.sdmnn/sig)then
         lbb1=lbnn1
         lbb2=lbnn2
         xmb1=xmnn1
         xmb2=xmnn2
      elseif(x1.le.(sdmnn+sdmnd)/sig)then
         lbb1=lbnd1
         lbb2=lbnd2
         xmb1=xmnd1
         xmb2=xmnd2
      elseif(x1.le.(sdmnn+sdmnd+sdmns)/sig)then
         lbb1=lbns1
         lbb2=lbns2
         xmb1=xmns1
         xmb2=xmns2
      elseif(x1.le.(sdmnn+sdmnd+sdmns+sdmnp)/sig)then
         lbb1=lbnp1
         lbb2=lbnp2
         xmb1=xmnp1
         xmb2=xmnp2
      elseif(x1.le.(sdmnn+sdmnd+sdmns+sdmnp+sdmdd)/sig)then
         lbb1=lbdd1
         lbb2=lbdd2
         xmb1=xmdd1
         xmb2=xmdd2
      elseif(x1.le.(sdmnn+sdmnd+sdmns+sdmnp+sdmdd+sdmds)/sig)then
         lbb1=lbds1
         lbb2=lbds2
         xmb1=xmds1
         xmb2=xmds2
      elseif(x1.le.(sdmnn+sdmnd+sdmns+sdmnp+sdmdd+sdmds+sdmdp)/sig)then
         lbb1=lbdp1
         lbb2=lbdp2
         xmb1=xmdp1
         xmb2=xmdp2
      elseif(x1.le.(sdmnn+sdmnd+sdmns+sdmnp+sdmdd+sdmds+sdmdp
     1        +sdmss)/sig)then
         lbb1=lbss1
         lbb2=lbss2
         xmb1=xmss1
         xmb2=xmss2
      elseif(x1.le.(sdmnn+sdmnd+sdmns+sdmnp+sdmdd+sdmds+sdmdp
     1        +sdmss+sdmsp)/sig)then
         lbb1=lbsp1
         lbb2=lbsp2
         xmb1=xmsp1
         xmb2=xmsp2
      elseif(x1.le.(sdmnn+sdmnd+sdmns+sdmnp+sdmdd+sdmds+sdmdp
     1        +sdmss+sdmsp+sdmpp)/sig)then
         lbb1=lbpp1
         lbb2=lbpp2
         xmb1=xmpp1
         xmb2=xmpp2
      else
         lbb1=lb1
         lbb2=lb2
         xmb1=em1
         xmb2=em2
         IBLOCK=504
      endif
      LB(I1)=lbb1
      E(i1)=xmb1
      LB(I2)=lbb2
      E(I2)=xmb2
      lb1=lb(i1)
      lb2=lb(i2)
      scheck=(s-(xmb1+xmb2)**2)*(s-(xmb1-xmb2)**2)
      if(scheck.lt.0) then
         write(99,*) 'scheck52: ', scheck
         scheck=0.
      endif
      pfinal=sqrt(scheck)/2./srt
      if(iblock.eq.502) then
         CALL dmangle(pxn,pyn,pzn,nt,ianti,pfinal,lbm)
      elseif(iblock.eq.504) then
         if(ianti.eq.0) then
            write (91,*) ' d+',lbm,' (regular d M elastic) @evt#',
     1           iaevt,' @nt=',nt,' lb1,2=',lb1,lb2
         else
            write (91,*) ' d+',lbm,' (regular dbar M elastic) @evt#',
     1           iaevt,' @nt=',nt,' lb1,2=',lb1,lb2
         endif
         CALL dmelangle(pxn,pyn,pzn,pfinal)
      else
         print *, 'Wrong iblock number in crdmbb()'
         stop
      endif
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
