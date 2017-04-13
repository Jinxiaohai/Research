        subroutine nnkaon(irun,iseed,ictrl,i1,i2,iblock,
     &                                   srt,pcx,pcy,pcz,nchrg)
      PARAMETER      (MAXSTR=150001,MAXR=1)
      PARAMETER      (AKA=0.498)
      COMMON   /AA/  R(3,MAXSTR)
      COMMON   /BB/  P(3,MAXSTR)
      COMMON   /CC/  E(MAXSTR)
      COMMON   /EE/  ID(MAXSTR),LB(MAXSTR)
      COMMON   /BG/BETAX,BETAY,BETAZ,GAMMA
      COMMON   /NN/NNN
      COMMON   /RUN/NUM
      COMMON   /PA/RPION(3,MAXSTR,MAXR)
      COMMON   /PB/PPION(3,MAXSTR,MAXR)
      COMMON   /PC/EPION(MAXSTR,MAXR)
      COMMON   /PD/LPION(MAXSTR,MAXR)
      dimension px(4),py(4),pz(4)
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
      SAVE   
      dm3=0.938
      dm4=0.938
      n=0
        if(nchrg.le.-1.or.nchrg.ge.3) dm3=1.232
        if(nchrg.eq.-2.or.nchrg.eq.4) dm4=1.232
          iblock = 0 
        call fstate(iseed,srt,dm3,dm4,px,py,pz,iflag)
        if(iflag.lt.0) then
           ictrl = -1
           n=n+1
           return
        endif
        iblock = 12
        pxrota=px(1)
        pyrota=py(1)
        pzrota=pz(1)
        call rotate(pcx,pcy,pcz,pxrota,pyrota,pzrota)
        px(1)=pxrota
        py(1)=pyrota
        pz(1)=pzrota
        pxrota=px(2)
        pyrota=py(2)
        pzrota=pz(2)
        call rotate(pcx,pcy,pcz,pxrota,pyrota,pzrota)
        px(2)=pxrota
        py(2)=pyrota
        pz(2)=pzrota
        pxrota=px(3)
        pyrota=py(3)
        pzrota=pz(3)
        call rotate(pcx,pcy,pcz,pxrota,pyrota,pzrota)
        px(3)=pxrota
        py(3)=pyrota
        pz(3)=pzrota
        pxrota=px(4)
        pyrota=py(4)
        pzrota=pz(4)
        call rotate(pcx,pcy,pcz,pxrota,pyrota,pzrota)
        px(4)=pxrota
        py(4)=pyrota
        pz(4)=pzrota
        nnn=nnn+2
        lpion(nnn,irun)=23
        if(nchrg.eq.-1.or.nchrg.eq.-2) then
        endif
        epion(nnn,irun)=aka
        lpion(nnn-1,irun)=21
        epion(nnn-1,irun)=aka
        e1cm   = sqrt(dm3**2 + px(1)**2 + py(1)**2 + pz(1)**2)
        p1beta = px(1)*betax + py(1)*betay + pz(1)*betaz
        transf = gamma * ( gamma*p1beta / (gamma+1) + e1cm)
        pt1i1 = betax*transf + px(1)
        pt2i1 = betay*transf + py(1)
        pt3i1 = betaz*transf + pz(1)
        eti1  = dm3
        lb1   = 2
        if(nchrg.ge.-2.and.nchrg.le.1) lb1=2
        if (nchrg .eq. -2 .or. nchrg .eq. -1) then
           lb1 = 6
        end if
        if(nchrg.eq.1.or.nchrg.eq.2) lb1=1
        if(nchrg.eq.3.or.nchrg.eq.4) lb1=9
        e2cm   = sqrt(dm4**2 + px(2)**2 + py(2)**2 + pz(2)**2)
        p2beta = px(2)*betax + py(2)*betay + pz(2)*betaz
        transf = gamma * ( gamma*p2beta / (gamma+1) + e2cm)
        pt1i2 = betax*transf + px(2)
        pt2i2 = betay*transf + py(2)
        pt3i2 = betaz*transf + pz(2)
        eti2  = dm4
        lb2   = 2
        if(nchrg.ge.-1.or.nchrg.le.1) lb2=2
        if(nchrg.eq. 2.or.nchrg.eq.3) lb2=1
        if(nchrg.eq. 4) lb2=9
        if(nchrg.eq.-2) lb2=6
                p(1,i1)=pt1i1
                p(2,i1)=pt2i1
                p(3,i1)=pt3i1
                e(i1)=eti1
                lb(i1)=lb1
                p(1,i2)=pt1i2
                p(2,i2)=pt2i2
                p(3,i2)=pt3i2
                e(i2)=eti2
                lb(i2)=lb2
        epcmk = sqrt(epion(nnn-1,irun)**2 + px(3)**2+py(3)**2+pz(3)**2)
        betak = px(3)*betax + py(3)*betay + pz(3)*betaz
        transf= gamma*(gamma*betak/(gamma+1.) + epcmk)
        ppion(1,nnn-1,irun)=betax*transf + px(3)
        ppion(2,nnn-1,irun)=betay*transf + py(3)
        ppion(3,nnn-1,irun)=betaz*transf + pz(3)
        rpion(1,nnn-1,irun)=r(1,i1)
        rpion(2,nnn-1,irun)=r(2,i1)
        rpion(3,nnn-1,irun)=r(3,i1)
        dppion(nnn-1,irun)=dpertp(i1)*dpertp(i2)
        epcmak = sqrt(epion(nnn,irun)**2 + px(4)**2 +py(4)**2+pz(4)**2)
        betaak = px(4)*betax + py(4)*betay + pz(4)*betaz
        transf= gamma*(gamma*betaak/(gamma+1.) + epcmak)
        ppion(1,nnn,irun)=betax*transf + px(4)
        ppion(2,nnn,irun)=betay*transf + py(4)
        ppion(3,nnn,irun)=betaz*transf + pz(4)
        rpion(1,nnn,irun)=r(1,i2)
        rpion(2,nnn,irun)=r(2,i2)
        rpion(3,nnn,irun)=r(3,i2)
        dppion(nnn,irun)=dpertp(i1)*dpertp(i2)
        return
        end
