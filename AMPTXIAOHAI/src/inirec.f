        subroutine inirec
        implicit double precision (a-h, o-z)
        external ran1
        parameter (MAXPTN=400001)
        common /para1/ mul
        common /para4/ iftflg, ireflg, igeflg, ibstfg
        common /para5/ iconfg, iordsc
        COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &       PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &       XMASS0(MAXPTN), ITYP0(MAXPTN)
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
        common /prec3/gxs(MAXPTN),gys(MAXPTN),gzs(MAXPTN),fts(MAXPTN),
     &     pxs(MAXPTN), pys(MAXPTN), pzs(MAXPTN), es(MAXPTN),
     &     xmasss(MAXPTN), ityps(MAXPTN)
        common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
        common /prec5/ eta(MAXPTN), rap(MAXPTN), tau(MAXPTN)
        common /prec6/ etas(MAXPTN), raps(MAXPTN), taus(MAXPTN)
        common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
        COMMON /ilist7/ LSTRG0(MAXPTN), LPART0(MAXPTN)
        COMMON /ilist8/ LSTRG1(MAXPTN), LPART1(MAXPTN)
        COMMON /smearz/smearp,smearh
        common /precpb/vxp(MAXPTN),vyp(MAXPTN),vzp(MAXPTN)
        common /precpa/vxp0(MAXPTN),vyp0(MAXPTN),vzp0(MAXPTN),
     1       xstrg0(MAXPTN),ystrg0(MAXPTN),
     2       xstrg(MAXPTN),ystrg(MAXPTN),istrg0(MAXPTN),istrg(MAXPTN)
        common/anim/nevent,isoft,isflag,izpc
        common /frzprc/ 
     &       gxfrz(MAXPTN), gyfrz(MAXPTN), gzfrz(MAXPTN), ftfrz(MAXPTN),
     &       pxfrz(MAXPTN), pyfrz(MAXPTN), pzfrz(MAXPTN), efrz(MAXPTN),
     &       xmfrz(MAXPTN), 
     &       tfrz(302), ifrz(MAXPTN), idfrz(MAXPTN), itlast
        common /rndm3/ iseedp
        common /para7/ ioscar,nsmbbbar,nsmmeson
        COMMON /AREVT/ IAEVT, IARUN, MISS
        SAVE   
        iseed=iseedp
        if(isoft.eq.5) then
           itlast=0
           call inifrz
        endif
        do 1001 i = 1, mul
           indxi=indx(i)
           ityp(i) = ityp0(indxi)
           gx(i) = gx0(indxi)
           gy(i) = gy0(indxi)
           gz(i) = gz0(indxi)
           ft(i) = ft0(indxi)
           px(i) = px0(indxi)
           py(i) = py0(indxi)
           pz(i) = pz0(indxi)
           e(i) = e0(indxi)
           xmass(i) = xmass0(indxi)
           LSTRG1(I) = LSTRG0(INDXI)
           LPART1(I) = LPART0(INDXI)
           vxp(I) = vxp0(INDXI)
           vyp(I) = vyp0(INDXI)
           vzp(I) = vzp0(INDXI)
           xstrg0(I) = xstrg(INDXI)
           ystrg0(I) = ystrg(INDXI)
           istrg0(I) = istrg(INDXI)
         if(isoft.eq.5) then
            idfrz(i)=ityp(i)
            gxfrz(i)=gx(i)
            gyfrz(i)=gy(i)
            gzfrz(i)=gz(i)
            ftfrz(i)=ft(i)
            pxfrz(i)=px(i)
            pyfrz(i)=py(i)
            pzfrz(i)=pz(i)
            efrz(i)=e(i)
            xmfrz(i)=xmass(i)
            ifrz(i)=0
         endif
 1001 continue
        do 1002 i = 1, mul
           ityps(i) = ityp(i)
           gxs(i) = gx(i)
           gys(i) = gy(i)
           gzs(i) = gz(i)
           fts(i) = ft(i)
           pxs(i) = px(i)
           pys(i) = py(i)
           pzs(i) = pz(i)
           es(i) = e(i)
           xmasss(i) = xmass(i)
 1002   continue
        if(isoft.eq.1.and.(ioscar.eq.2.or.ioscar.eq.3))
     1       write(92,*) iaevt,miss,mul
        do 1003 i = 1, mul
           energy = e(i)
           vx(i) = px(i) / energy
           vy(i) = py(i) / energy
           vz(i) = pz(i) / energy
           if (iftflg .eq. 0) then
              formt = ft(i)
            if(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) then
               gx(i) = gx(i) + vxp(i) * formt
               gy(i) = gy(i) + vyp(i) * formt
               gz(i) = gz(i) + vzp(i) * formt
            else
               gx(i) = gx(i) + vx(i) * formt
               gy(i) = gy(i) + vy(i) * formt
               gz(i) = gz(i) + vz(i) * formt
            endif
           end if
           if(ioscar.eq.2.or.ioscar.eq.3) then
              if(dmax1(abs(gx(i)),abs(gy(i)),
     1             abs(gz(i)),abs(ft(i))).lt.9999) then
                 write(92,200) ityp(i),px(i),py(i),pz(i),xmass(i),
     1           gx(i),gy(i),gz(i),ft(i),istrg0(i),xstrg0(i),ystrg0(i)
              else
                 write(92,201) ityp(i),px(i),py(i),pz(i),xmass(i),
     1           gx(i),gy(i),gz(i),ft(i),istrg0(i),xstrg0(i),ystrg0(i)
              endif
           endif
 200       format(I3,2(1x,f7.2),1x,f8.2,1x,f6.3,4(1x,f8.2),
     1          1x,I5,2(1x,f7.2))
 201       format(I3,2(1x,f7.2),1x,f8.2,1x,f6.3,4(1x,e8.2),
     1          1x,I5,2(1x,f7.2))
 1003   continue
        if (iconfg .le. 3) then
           do 1004 i = 1, mul
              if (ft(i) .le. abs(gz(i))) then
                 eta(i) = 1000000.d0
              else
                 eta(i) = 0.5d0 * log((ft(i) + gz(i)) / (ft(i) - gz(i)))
              end if
              if (e(i) .le. abs(pz(i))) then
                 rap(i) = 1000000.d0
              else
                 rap(i) = 0.5d0 * log((e(i) + pz(i)) / (e(i) - pz(i)))
              end if
              if(eta(i).lt.1000000.d0) then
                 tau(i) = ft(i) / cosh(eta(i))
              else
                 tau(i) = 1d-10
              endif
 1004      continue
           do 1005 i = 1, mul
              etas(i) = eta(i)
              raps(i) = rap(i)
              taus(i) = tau(i)
 1005      continue
        end if
        return
        end
