        subroutine newka(icase,irun,iseed,dt,nt,ictrl,i1,i2,
     &                                   srt,pcx,pcy,pcz,iblock)
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
      COMMON/RNDF77/NSEED
      SAVE   
        logical lb1bn, lb2bn,lb1mn,lb2mn
        logical lb1bn1, lb2bn1, lb1bn0, lb2bn0
        logical lb1mn0, lb2mn0, lb1mn1, lb2mn1
        logical lb1mn2, lb2mn2
        icase=-1
        nchrg=-100
        ictrl = 1
        lb1=lb(i1)
        lb2=lb(i2)
        em1=e(i1)
        em2=e(i2)
        lb1bn=lb1.eq.1.or.lb1.eq.2.or.(lb1.gt.5.and.lb1.le.13)
        lb2bn=lb2.eq.1.or.lb2.eq.2.or.(lb2.gt.5.and.lb2.le.13)
        lb1bn0=lb1.eq.2.or.lb1.eq.7.or.lb1.eq.10.or.lb1.eq.12
        lb2bn0=lb2.eq.2.or.lb2.eq.7.or.lb2.eq.10.or.lb2.eq.12
        lb1bn1=lb1.eq.1.or.lb1.eq.8.or.lb1.eq.11.or.lb1.eq.13
        lb2bn1=lb2.eq.1.or.lb2.eq.8.or.lb2.eq.11.or.lb2.eq.13
        lb1mn=em1.lt.0.2.or.lb1.eq.0.or.(lb1.ge.25.and.lb1.le.29)
        lb2mn=em2.lt.0.2.or.lb2.eq.0.or.(lb2.ge.25.and.lb2.le.29)
        lb1mn0=lb1.eq.0.or.lb1.eq.4.or.lb1.eq.26.or.
     &                        lb1.eq.28.or.lb1.eq.29
        lb2mn0=lb2.eq.0.or.lb2.eq.4.or.lb2.eq.26.or.
     &                        lb2.eq.28.or.lb2.eq.29
        lb1mn1= lb1.eq.5.or.lb1.eq.27
        lb2mn1= lb2.eq.5.or.lb2.eq.27
        lb1mn2=lb1.eq.3.or.lb1.eq.25
        lb2mn2=lb2.eq.3.or.lb2.eq.25
        if(lb1bn.and.lb2bn) then
           icase=1
           sig=40.
           if(lb1.eq.9.and.lb2.eq.9) then
                nchrg=4
           endif   
           if((lb1bn1.and.lb2.eq.9)
     &        .or.(lb2bn1.and.lb1.eq.9))then
                nchrg=3
           endif
           if((lb1bn0.and.lb2.eq.9)
     &        .or.(lb2bn0.and.lb1.eq.9)
     &        .or.(lb1bn1.and.lb2bn1)) then
                   nchrg=2
           endif
           if((lb1bn1.and.lb2bn0).or.(lb1.eq.6.and.lb2.eq.9)
     &        .or.(lb2bn1.and.lb1bn0)
     &        .or.(lb2.eq.6.and.lb1.eq.9))then
                   nchrg=1
           endif
           if((lb1bn0.and.lb2bn0).or.(lb1bn1.and.lb2.eq.6)
     &              .or.(lb2bn1.and.lb1.eq.6)) then
                   nchrg=0
           endif
           if((lb1bn0.and.lb2.eq.6)
     &        .or.(lb2bn0.and.lb1.eq.6))then
                nchrg=-1
           endif
           if(lb1.eq.6.and.lb2.eq.6) then
                nchrg=-2
           endif
           if(nchrg.ge.-1.and.nchrg.le.2) then
                   brsig = x2kaon(srt)
           else
                   brsig=0.0
           endif
           BRSIG = 2.0 * BRSIG
        endif
        if((lb1bn.and.lb2mn).OR.(lb2bn.and.lb1mn)) then
          icase=2
          sig=20.
          sigma0 = piNsg0(srt)
          brsig=0.0
          if((lb1bn1.and.lb2mn0)
     &       .or.(lb2bn1.and.lb1mn0).
     & or.(lb1bn0.and.lb2mn1).or.(lb2bn0.and.lb1mn1).
     & or.(lb1.eq.9.and.lb2mn2).or.(lb2.eq.9.and.lb1mn2))then
                nchrg=1
                if(lb1bn1.or.lb2bn1) brsig=0.5*sigma0
                if(lb1bn0.or.lb2bn0) brsig=2.0*sigma0
          endif
          if( (lb1bn0.and.lb2mn0 )
     &       .or.(lb2bn0.and.lb1mn0)
     &  .or.(lb1bn1.and.lb2mn2).or.(lb2bn1.and.lb1mn2)
     &  .or.(lb1.eq.6.and.lb2mn1).or.(lb2.eq.6.and.lb1mn1)) then
                nchrg=0
                if(lb1bn1.or.lb2bn1) then
                  brsig=3.0*sigma0
                  ratiok = 2./3.
                endif
                if(lb1bn0.or.lb2bn0) then
                  brsig=2.5*sigma0
                  ratiok = 0.2
                endif
          endif
          if( (lb1bn0.and.lb2mn2)
     &       .or.(lb2bn0.and.lb1mn2)
     & .or.(lb1.eq.6.and.lb2mn0).or.(lb2.eq.6.and.lb1mn0)) then
                nchrg=-1
                if(lb1bn0.or.lb2bn0) brsig=sigma0
          endif
          if((lb1.eq.6.and.lb2mn2)
     &       .or.(lb2.eq.6.and.lb1mn2))then
                nchrg=-2
          endif
          if((lb1bn1.and.lb2mn1)
     &       .or.(lb2bn1.and.lb1mn1)
     & .or.(lb1.eq.9.and.lb2mn0).or.(lb2.eq.9.and.lb1mn0)) then
                nchrg=2
          endif
          IF (NCHRG .GE. -2 .AND. NCHRG .LE. 2) THEN
             BRSIG = 3.0 * SIGMA0
          END IF
        endif
        if( (lb1bn.and.(lb2.eq.21.or.lb2.eq.-30)).OR.
     &     (lb2bn.and.(lb1.eq.21.or.lb1.eq.-30)) )then 
          bmass=0.938
          if(srt.le.(bmass+aka)) then
                pkaon=0.
          else
            pkaon=sqrt(((srt**2-(aka**2+bmass**2))/2./bmass)**2-aka**2)
          endif
          sig=0.
          if(lb1.eq.1.or.lb2.eq.1.or.lb1.eq.8.or.lb2.eq.8.or.
     &    lb1.eq.11.or.lb2.eq.11.or.lb1.eq.13.or.lb2.eq.13) then
              nchrg=0
              sigela=akPel(pkaon)
              sigsgm=3.*akPsgm(pkaon)
              sig=sigela+sigsgm+akPlam(pkaon)
          endif
          if(lb1.eq.2.or.lb2.eq.2.or.lb1.eq.7.or.lb2.eq.7.or.
     &    lb1.eq.10.or.lb2.eq.10.or.lb1.eq.12.or.lb2.eq.12) then
              nchrg=-1
              sigela=akNel(pkaon)
              sigsgm=2.*akNsgm(pkaon)
              sig=sigela+sigsgm+akNlam(pkaon)
          endif
          if(lb1.eq.6.or.lb2.eq.6) then
              nchrg=-2
              sigela=akNel(pkaon)
              sigsgm=akNsgm(pkaon)
              sig=sigela+sigsgm
          endif
          if(lb1.eq.9.or.lb2.eq.9) then
              nchrg=1
              sigela=akPel(pkaon)
              sigsgm=2.*akPsgm(pkaon)
              sig=sigela+sigsgm+akPlam(pkaon)
          endif
          sigela = 0.5 * (AKPEL(PKAON) + AKNEL(PKAON))
          SIGSGM = 1.5 * AKPSGM(PKAON) + AKNSGM(PKAON)
          SIG = sigela + SIGSGM + AKPLAM(PKAON)
          if(sig.gt.1.e-7) then
              icase=3
              brel=sigela/sig
              brsgm=sigsgm/sig
              brsig = sig
          endif
        endif
        if(((lb1.ge.14.and.lb1.le.17).and.(lb2.ge.3.and.lb2.le.5)).OR.
     &     ((lb2.ge.14.and.lb2.le.17).and.(lb1.ge.3.and.lb1.le.5)))then
           nchrg=-100
           if((lb1.eq.15.and.(lb2.eq.3.or.lb2.eq.25)).OR.
     &              (lb2.eq.15.and.(lb1.eq.3.or.lb1.eq.25))) then
                nchrg=-2
                  bmass=1.232
           endif
           if((lb1.eq.15.and.lb2mn0).or.(lb2.eq.15.and.lb1mn0).OR.
     &       ((lb1.eq.14.or.lb1.eq.16).and.(lb2.eq.3.or.lb2.eq.25)).OR.
     &       ((lb2.eq.14.or.lb2.eq.16).and.(lb1.eq.3.or.lb1.eq.25)))then
                nchrg=-1
                 bmass=0.938
           endif
           if((lb1.eq.15.and.(lb2.eq.5.or.lb2.eq.27)).OR.
     &              (lb2.eq.15.and.(lb1.eq.5.or.lb1.eq.27)).or.
     &        (lb1.eq.17.and.(lb2.eq.3.or.lb2.eq.25)).OR.
     &              (lb2.eq.17.and.(lb1.eq.3.or.lb1.eq.25)).or.
     &       ((lb1.eq.14.or.lb1.eq.16).and.lb2mn0).OR.
     &       ((lb2.eq.14.or.lb2.eq.16).and.lb1mn0)) then
                nchrg=0
                 bmass=0.938
           endif
           if((lb1.eq.17.and.lb2mn0).or.(lb2.eq.17.and.lb1mn0).OR.
     &       ((lb1.eq.14.or.lb1.eq.16).and.(lb2.eq.5.or.lb2.eq.27)).OR.
     &       ((lb2.eq.14.or.lb2.eq.16).and.(lb1.eq.5.or.lb1.eq.27)))then
                nchrg=1
                 bmass=1.232
           endif
           sig = 0.
           if(nchrg.ne.-100.and.srt.gt.(aka+bmass)) then
             icase=4
             pkaon=sqrt(((srt**2-(aka**2+0.938**2))/2./0.938)**2-aka**2)
             if(lb1.eq.14.or.lb2.eq.14) then
                if(nchrg.ge.0) sigma0=akPlam(pkaon)
                if(nchrg.lt.0) sigma0=akNlam(pkaon)
             else
                if(nchrg.ge.0) sigma0=akPsgm(pkaon)
                if(nchrg.lt.0) sigma0=akNsgm(pkaon)
                SIGMA0 = 1.5 * AKPSGM(PKAON) + AKNSGM(PKAON)
             endif
             sig=(srt**2-(aka+bmass)**2)*(srt**2-(aka-bmass)**2)/
     &         (srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)*sigma0
             if(nchrg.eq.-2.or.nchrg.eq.2) sig=2.*sig
             IF (LB1 .EQ. 14 .OR. LB2 .EQ. 14) THEN
                SIG = 4.0 / 3.0 * SIG
             ELSE IF (NCHRG .EQ. -2 .OR. NCHRG .EQ. 2) THEN
                SIG = 8.0 / 9.0 * SIG
             ELSE
                SIG = 4.0 / 9.0 * SIG
             END IF
             brsig = sig
             if(sig.lt.1.e-7) sig = 1.e-7
           endif
           icase=4
           brsig = sig
           sigela = 10.
           sig = sig + sigela
           brel = sigela/sig
        endif
        if(icase.eq.-1) then
           ictrl = -1
           return
        endif
        px1cm=pcx
        py1cm=pcy
        pz1cm=pcz
        ds=sqrt(sig/31.4)
        dsr=ds+0.1
        ec=(em1+em2+0.02)**2
        call distce(i1,i2,dsr,ds,dt,ec,srt,ic,px1cm,py1cm,pz1cm)
        if(ic.eq.-1) then
           ictrl = -1
           return
        endif
        ik=0
        ik0=0
        ik1=0
        ik2=0
        ik3=0
        il=0
        im=0
        im3=0
        im4=0
        in=0
        inpion=0
        ipipi=0
        sgsum=0.
        sgsum1=0.
        sgsum3=0.
        if(icase.eq.1) then
           ik=ik+1
           if(srt.gt.2.8639) then
                ik0=ik0+1
                if(em1.lt.1.0.and.em2.lt.1.0) then
                        ik1=ik1+1
                        sgsum1=sgsum1+brsig
                endif
                if(em1.gt.1.0.and.em2.gt.1.0) then
                        ik3=ik3+1
                        sgsum3=sgsum3+brsig
                endif
                if(em1.gt.1.0.and.em2.lt.1.0) ik2=ik2+1
                if(em1.lt.1.0.and.em2.gt.1.0) ik2=ik2+1
                sgsum=sgsum+brsig
           endif
        endif
        if(icase.eq.2) inpion=inpion+1
        if(icase.eq.5) ipipi=ipipi+1
        if(RANART(NSEED).gt.(brsig/sig)) then
           ictrl = -1
           return
        endif
        il=il+1
        if(icase.eq.1) then
          in=in+1
          call nnkaon(irun,iseed,
     &          ictrl,i1,i2,iblock,srt,pcx,pcy,pcz,nchrg)
        endif
        if(icase.eq.2) then
          im=im+1
          call npik(irun,iseed,dt,nt,ictrl,i1,i2,srt,
     &              pcx,pcy,pcz,nchrg,ratiok,iblock)
        endif
        if(icase.eq.3) then
          im3=im3+1
          call kaonN(brel,brsgm,irun,iseed,dt,nt,ictrl,
     &                i1,i2,iblock,srt,pcx,pcy,pcz,nchrg)
        endif
        if(icase.eq.4) then
          im4=im4+1
      if(RANART(NSEED).lt.brel) then
         ielstc = 1
      else
         ielstc = 0
      endif                  
          call Pihypn(ielstc,irun,iseed,dt,nt,ictrl,i1,i2,srt,
     &                   pcx,pcy,pcz,nchrg,iblock)
        endif
        return
        end
