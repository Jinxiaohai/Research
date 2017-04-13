        subroutine zpcrun(*)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        parameter (tend1 = 250d0)
        parameter (tend2 = 6.1d0)
        common /para1/ mul
        common /para5/ iconfg, iordsc
        common /para7/ ioscar,nsmbbbar,nsmmeson
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
        common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
        common /prec5/ eta(MAXPTN), rap(MAXPTN), tau(MAXPTN)
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
        common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
        common /ilist6/ t, iopern, icolln
        common /ana1/ ts(12)
        common/anim/nevent,isoft,isflag,izpc
        COMMON /AREVT/ IAEVT, IARUN, MISS
        SAVE   
        if (mod(ictype, 2) .eq. 0) then
           call savrec(iscat)
           call savrec(jscat)
        end if
        call getict(t1)
        if (iconfg .eq. 1 .and. t1 .gt. tlarge / 2d0) return 1
        if (iconfg .eq. 2 .or. iconfg .eq. 3) then
           if (t1 .gt. tend1) return 1
        end if
        if (iconfg .eq. 4 .or. iconfg .eq. 5) then
           if (t1 .gt. tend2) return 1
        end if
      if(isoft.eq.5) then
         call local(t1)
      endif
        iopern = iopern + 1
        t = t1
        if (mod(ictype, 2) .eq. 0) then
           icolln = icolln + 1
        end if
        if (iconfg .eq. 1
     &     .or. iconfg .eq. 2
     &     .or. iconfg .eq. 4) then
           if (ictype .eq. 1 .or. ictype .eq. 2 .or. 
     &        ictype .eq. 5 .or. ictype .eq. 6) then
              call celasn
           end if
        end if
        if (ictype .ne. 1) then
           iscat0 = iscat
           jscat0 = jscat
           iscat = max0(iscat0, jscat0)
           jscat = min0(iscat0, jscat0)
           if (jscat .ne. 0) then
              if(next(jscat) .ne. iscat) then
                 print *, 'iscat=', iscat, 'jscat=', jscat,
     &                'next(', jscat, ')=', next(jscat)
                 if (ct(iscat) .lt. tlarge / 2d0) stop 'tterr'
                 if (ct(jscat) .lt. tlarge / 2d0) stop 'tterr'
              endif
           end if 
         niscat=iscat
         njscat=jscat
           if (icsta(iscat) .ne. 0) call cellre(niscat, t)
           if (jscat .ne. 0) then
              if (icsta(jscat) .ne. 0) call cellre(njscat, t)
           end if
           if (mod(ictype, 2) .eq. 0) then
              if(ioscar.eq.3) then
            write(95,*) 'event,miss,iscat,jscat=',iaevt,miss,iscat,jscat
                 if(dmax1(abs(gx(iscat)),abs(gy(iscat)),
     1                abs(gz(iscat)),abs(ft(iscat)),abs(gx(jscat)),
     2                abs(gy(jscat)),abs(gz(jscat)),abs(ft(jscat)))
     3                .lt.9999) then
                    write(95,200) ityp(iscat),px(iscat),py(iscat),
     1                   pz(iscat),xmass(iscat),gx(iscat),gy(iscat),
     2                   gz(iscat),ft(iscat)
                    write(95,200) ityp(jscat),px(jscat),py(jscat),
     1                   pz(jscat),xmass(jscat),gx(jscat),gy(jscat),
     2                   gz(jscat),ft(jscat)
                 else
                    write(95,201) ityp(iscat),px(iscat),py(iscat),
     1                   pz(iscat),xmass(iscat),gx(iscat),gy(iscat),
     2                   gz(iscat),ft(iscat)
                    write(95,201) ityp(jscat),px(jscat),py(jscat),
     1                   pz(jscat),xmass(jscat),gx(jscat),gy(jscat),
     2                   gz(jscat),ft(jscat)
                 endif
              endif
              call scat(t, iscat, jscat)
              if(ioscar.eq.3) then
                 if(dmax1(abs(gx(iscat)),abs(gy(iscat)),
     1                abs(gz(iscat)),abs(ft(iscat)),abs(gx(jscat)),
     2                abs(gy(jscat)),abs(gz(jscat)),abs(ft(jscat)))
     3                .lt.9999) then
                    write(95,200) ityp(iscat),px(iscat),py(iscat),
     1                   pz(iscat),xmass(iscat),gx(iscat),gy(iscat),
     2                   gz(iscat),ft(iscat)
                    write(95,200) ityp(jscat),px(jscat),py(jscat),
     1                   pz(jscat),xmass(jscat),gx(jscat),gy(jscat),
     2                   gz(jscat),ft(jscat)
                 else
                    write(95,201) ityp(iscat),px(iscat),py(iscat),
     1                   pz(iscat),xmass(iscat),gx(iscat),gy(iscat),
     2                   gz(iscat),ft(iscat)
                    write(95,201) ityp(jscat),px(jscat),py(jscat),
     1                   pz(jscat),xmass(jscat),gx(jscat),gy(jscat),
     2                   gz(jscat),ft(jscat)
                 endif
              endif
           endif
 200       format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f8.2))
 201       format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e8.2))
        end if
        call ulist(t)
        if (ifmpt .le. mul) then
           if (ictype .ne. 0 .and. ictype .ne. 3 
     &        .and. ictype .ne. 4) then
              ichkpt = ichkpt + 1
              ifmpt = ifmpt + 1
           end if
        end if
        return
        end
