        subroutine flowp(idd)
        implicit double precision  (a-h, o-z)
        real dt
        parameter (MAXPTN=400001)
        parameter (bmt=0.05d0)
        dimension nlfile(3),nsfile(3),nmfile(3)
        dimension v2pp(3),xnpp(3),v2psum(3),v2p2sm(3),nfile(3)
        dimension tsp(31),v2pevt(3),v2pavg(3),varv2p(3)
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
        COMMON /para1/ mul
        COMMON /prec2/GX5(MAXPTN),GY5(MAXPTN),GZ5(MAXPTN),FT5(MAXPTN),
     &       PX5(MAXPTN), PY5(MAXPTN), PZ5(MAXPTN), E5(MAXPTN),
     &       XMASS5(MAXPTN), ITYP5(MAXPTN)
        COMMON /pflow/ v2p(30,3),xnpart(30,3),etp(30,3),
     1       s2p(30,3),v2p2(30,3),nevt(30)
        COMMON /pflowf/ v2pf(30,3),xnpf(30,3),etpf(30,3),
     1                 xncoll(30),s2pf(30,3),v2pf2(30,3)
        COMMON /pfrz/ v2pfrz(30,3),xnpfrz(30,3),etpfrz(30,3),
     1       s2pfrz(30,3),v2p2fz(30,3),tscatt(31),
     2       nevtfz(30),iscatt(30)
        COMMON /hflow/ v2h(30,3),xnhadr(30,3),eth(30,3),
     1 v2h2(30,3),s2h(30,3)
        COMMON /AREVT/ IAEVT, IARUN, MISS
        common/anim/nevent,isoft,isflag,izpc
        common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
        COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &   IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
        common /precpb/vxp(MAXPTN),vyp(MAXPTN),vzp(MAXPTN)
        SAVE   
        dimension etpl(30,3),etps(30,3),etplf(30,3),etpsf(30,3),
     &       etlfrz(30,3),etsfrz(30,3),
     &       xnpl(30,3),xnps(30,3),xnplf(30,3),xnpsf(30,3),
     &       xnlfrz(30,3),xnsfrz(30,3),
     &       v2pl(30,3),v2ps(30,3),v2plf(30,3),v2psf(30,3),
     &       s2pl(30,3),s2ps(30,3),s2plf(30,3),s2psf(30,3),
     &       DMYil(50,3),DMYfl(50,3),
     &       DMYis(50,3),DMYfs(50,3)
        data tsp/0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
     &       1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
     &       2  , 3,   4,   5,   6,   7,   8,   9,   10,  20,  30/
        if(idd.eq.0) then        
           nfile(1)=60
           nfile(2)=64
           nfile(3)=20
           OPEN (nfile(1),FILE='ana1/v2p.dat', STATUS = 'UNKNOWN')
           OPEN (nfile(1)+1, 
     1 FILE = 'ana1/v2p-formed.dat', STATUS = 'UNKNOWN')
           OPEN (nfile(1)+2, 
     1 FILE = 'ana1/v2p-active.dat', STATUS = 'UNKNOWN')
           OPEN (nfile(1)+3, 
     1 FILE = 'ana1/v2ph.dat', STATUS = 'UNKNOWN')
           OPEN (nfile(2),FILE='ana1/v2p-y2.dat', STATUS = 'UNKNOWN')
           OPEN (nfile(2)+1, 
     1 FILE = 'ana1/v2p-formed2.dat', STATUS = 'UNKNOWN')
           OPEN (nfile(2)+2, 
     1 FILE = 'ana1/v2p-active2.dat', STATUS = 'UNKNOWN')
           OPEN (nfile(2)+3, 
     1 FILE = 'ana1/v2ph-y2.dat', STATUS = 'UNKNOWN')
           OPEN (nfile(3),FILE='ana1/v2p-y1.dat', STATUS = 'UNKNOWN')
           OPEN (nfile(3)+1, 
     1 FILE = 'ana1/v2p-formed1.dat', STATUS = 'UNKNOWN')
           OPEN (nfile(3)+2, 
     1 FILE = 'ana1/v2p-active1.dat', STATUS = 'UNKNOWN')
           OPEN (nfile(3)+3, 
     1 FILE = 'ana1/v2ph-y1.dat', STATUS = 'UNKNOWN')
           OPEN (49, FILE = 'ana1/v2p-ebe.dat', STATUS = 'UNKNOWN')
           write(49, *) '    ievt,  v2p,  v2p_y2,   v2p_y1'
           OPEN (59, FILE = 'ana1/v2h.dat', STATUS = 'UNKNOWN')
           OPEN (68, FILE = 'ana1/v2h-y2.dat', STATUS = 'UNKNOWN')
           OPEN (69, FILE = 'ana1/v2h-y1.dat', STATUS = 'UNKNOWN')
           OPEN (88, FILE = 'ana1/v2h-ebe.dat', STATUS = 'UNKNOWN')
           write(88, *) '    ievt,  v2h,  v2h_y2,   v2h_y1'
           nlfile(1)=70
           nlfile(2)=72
           nlfile(3)=74
           OPEN (nlfile(1),FILE='ana1/mtl.dat', STATUS = 'UNKNOWN')
           OPEN (nlfile(1)+1, 
     1 FILE = 'ana1/mtl-formed.dat', STATUS = 'UNKNOWN')
           OPEN (nlfile(2),FILE='ana1/mtl-y2.dat', STATUS = 'UNKNOWN')
           OPEN (nlfile(2)+1, 
     1 FILE = 'ana1/mtl-formed2.dat', STATUS = 'UNKNOWN')
           OPEN (nlfile(3),FILE='ana1/mtl-y1.dat', STATUS = 'UNKNOWN')
           OPEN (nlfile(3)+1, 
     1 FILE = 'ana1/mtl-formed1.dat', STATUS = 'UNKNOWN')
           nsfile(1)=76
           nsfile(2)=78
           nsfile(3)=80
           OPEN (nsfile(1),FILE='ana1/mts.dat', STATUS = 'UNKNOWN')
           OPEN (nsfile(1)+1, 
     1 FILE = 'ana1/mts-formed.dat', STATUS = 'UNKNOWN')
           OPEN (nsfile(2),FILE='ana1/mts-y2.dat', STATUS = 'UNKNOWN')
           OPEN (nsfile(2)+1, 
     1 FILE = 'ana1/mts-formed2.dat', STATUS = 'UNKNOWN')
           OPEN (nsfile(3),FILE='ana1/mts-y1.dat', STATUS = 'UNKNOWN')
           OPEN (nsfile(3)+1, 
     1 FILE = 'ana1/mts-formed1.dat', STATUS = 'UNKNOWN')
           nmfile(1)=82
           nmfile(2)=83
           nmfile(3)=84
           OPEN (nmfile(1),FILE='ana1/Nmt.dat', STATUS = 'UNKNOWN')
           OPEN (nmfile(2),FILE='ana1/Nmt-y2.dat', STATUS = 'UNKNOWN')
           OPEN (nmfile(3),FILE='ana1/Nmt-y1.dat', STATUS = 'UNKNOWN')
           ifanim=0
           if(ifanim.eq.1) then
              OPEN (10, FILE = 'ana1/h-animate.dat', STATUS = 'UNKNOWN')
              write(10,*) ntmax, dt
              OPEN (11, FILE = 'ana1/p-animate.dat', STATUS = 'UNKNOWN')
              OPEN (15, FILE = 'ana1/p-finalft.dat', STATUS = 'UNKNOWN')
           endif
           if(nevent.lt.1) 
     1          OPEN (93, FILE = 'ana1/parton-t.dat', STATUS='UNKNOWN')
           itimep=0
           itanim=0
           iaevtp=0
           do 1002 ii=1,50
              do 1001 iy=1,3
                 DMYil(ii,iy) = 0d0
                 DMYfl(ii,iy) = 0d0
                 DMYis(ii,iy) = 0d0
                 DMYfs(ii,iy) = 0d0
 1001         continue
 1002      continue
           do 1003 ii=1,31
              tscatt(ii)=0d0
 1003      continue
           do 1005 ii=1,30
              nevt(ii)=0
              xncoll(ii)=0d0
              nevtfz(ii)=0
              iscatt(ii)=0
              do 1004 iy=1,3
                 v2p(ii,iy)=0d0
                 v2p2(ii,iy)=0d0
                 s2p(ii,iy)=0d0
                 etp(ii,iy)=0d0
                 xnpart(ii,iy)=0d0
                 v2pf(ii,iy)=0d0
                 v2pf2(ii,iy)=0d0
                 s2pf(ii,iy)=0d0
                 etpf(ii,iy)=0d0
                 xnpf(ii,iy)=0d0
                 v2pfrz(ii,iy)=0d0
                 v2p2fz(ii,iy)=0d0
                 s2pfrz(ii,iy)=0d0
                 etpfrz(ii,iy)=0d0
                 xnpfrz(ii,iy)=0d0
                 etpl(ii,iy)=0d0
                 etps(ii,iy)=0d0
                 etplf(ii,iy)=0d0
                 etpsf(ii,iy)=0d0
                 etlfrz(ii,iy)=0d0
                 etsfrz(ii,iy)=0d0
              xnpl(ii,iy)=0d0
              xnps(ii,iy)=0d0
              xnplf(ii,iy)=0d0
              xnpsf(ii,iy)=0d0
              xnlfrz(ii,iy)=0d0
              xnsfrz(ii,iy)=0d0
              v2pl(ii,iy)=0d0
              v2ps(ii,iy)=0d0
              v2plf(ii,iy)=0d0
              v2psf(ii,iy)=0d0
              s2pl(ii,iy)=0d0
              s2ps(ii,iy)=0d0
              s2plf(ii,iy)=0d0
              s2psf(ii,iy)=0d0
 1004      continue
 1005   continue
           do 1006 iy=1,3
              v2pevt(iy)=0d0
              v2pavg(iy)=0d0
              varv2p(iy)=0d0
              v2pp(iy)=0.d0
              xnpp(iy)=0d0
              v2psum(iy)=0.d0
              v2p2sm(iy)=0.d0
 1006      continue
        else if(idd.eq.1) then        
           if(iaevt.ne.iaevtp.and.ianp.eq.31) itanim=0
           t2time = FT5(iscat)
           do 1008 ianp = 1, 30
              if (t2time.lt.tsp(ianp+1).and.t2time.ge.tsp(ianp)) then
                 xncoll(ianp)=xncoll(ianp)+1d0
                 if(ianp.le.itimep.and.iaevt.eq.iaevtp) goto 101
                 nevt(ianp)=nevt(ianp)+1
                 tscatt(ianp+1)=t2time
                 iscatt(ianp)=1
                 nevtfz(ianp)=nevtfz(ianp)+1
                 do 100 I=1,mul
                    delta=1d-8
                    if((E5(i)-dabs(PZ5(i))+delta).le.0) then
                       ypartn=1000000.d0*sign(1.d0,PZ5(i))
                       write(6,*) 'ypartn error',E5(i)-dabs(PZ5(i))
                    else
                       ypartn=0.5d0*dlog((E5(i)+PZ5(i)+delta)
     1                      /(E5(i)-PZ5(i)+delta))
                    endif
                    pt2=PX5(I)**2+PY5(I)**2
                    iloop=1
                    if(dabs(ypartn).le.1d0) then
                       iloop=2
                       if(dabs(ypartn).le.0.5d0) then
                          iloop=3
                       endif
                    endif
                    do 50 iy=1,iloop
                       if(pt2.gt.0d0) then
                          v2prtn=(PX5(I)**2-PY5(I)**2)/pt2
                          if(dabs(v2prtn).gt.1d0) 
     1 write(nfile(iy),*) 'v2prtn>1',v2prtn
                          v2p(ianp,iy)=v2p(ianp,iy)+v2prtn
                          v2p2(ianp,iy)=v2p2(ianp,iy)+v2prtn**2
                       endif
                       xperp2=GX5(I)**2+GY5(I)**2
                       if(xperp2.gt.0d0) 
     1        s2p(ianp,iy)=s2p(ianp,iy)+(GX5(I)**2-GY5(I)**2)/xperp2
                       xnpart(ianp,iy)=xnpart(ianp,iy)+1d0
                       etp(ianp,iy)=etp(ianp,iy)+dsqrt(pt2+XMASS5(I)**2)
                       if(FT5(I).le.t2time) then
                          v2pf(ianp,iy)=v2pf(ianp,iy)+v2prtn
                          v2pf2(ianp,iy)=v2pf2(ianp,iy)+v2prtn**2
                          if(xperp2.gt.0d0) 
     1        s2pf(ianp,iy)=s2pf(ianp,iy)+(GX5(I)**2-GY5(I)**2)/xperp2
                          xnpf(ianp,iy)=xnpf(ianp,iy)+1d0
                  etpf(ianp,iy)=etpf(ianp,iy)+dsqrt(pt2+XMASS5(I)**2)
                       endif
 50                    continue
 100                 continue
                 itimep=ianp
                 iaevtp=iaevt
                 if(ianp.eq.30) then
                    do 1007 iy=1,3
                       npartn=IDINT(xnpart(ianp,iy)-xnpp(iy))
                       if(npartn.ne.0) then
                          v2pevt(iy)=(v2p(ianp,iy)-v2pp(iy))/npartn
                          v2psum(iy)=v2psum(iy)+v2pevt(iy)
                          v2p2sm(iy)=v2p2sm(iy)+v2pevt(iy)**2
                          v2pp(iy)=v2p(ianp,iy)
                          xnpp(iy)=xnpart(ianp,iy)
                       endif
 1007               continue
                    write(49, 160) iaevt,v2pevt
                 endif
                 goto 101
              endif
 1008      continue
 101       if(nevent.lt.1) then
              do 110 nt = 1, ntmax
                 time1=dble(nt*dt)
                 time2=dble((nt+1)*dt)
                 if (t2time.lt.time2.and.t2time.ge.time1) then
                    if(nt.le.itanim) return
                    if(ifanim.eq.1) write(11,*) t2time
                    iform=0
                    ne1all=0
                    ne1form=0
                    do 1009 I=1,mul
                       gz_now=GZ5(i)+(t2time-FT5(i))*PZ5(i)/E5(i)
                       If(dabs(gz_now).lt.t2time) then
                       etap=0.5d0*dlog((t2time+gz_now)/(t2time-gz_now))
                       else
                          etap=1000000.d0*sign(1.d0,gz_now)
                       endif
                       ne1all=ne1all+1
                       if(FT5(I).le.t2time) ne1form=ne1form+1
                       if(FT5(I).le.t2time) iform=iform+1
 1009               continue
                    if(ifanim.eq.1) write(11,*) iform
                    write(93,184) 'evt#,t,np,npformed=',
     1                   iaevt,t2time,ne1all,ne1form
 184                format(a20,i7,f8.4,2(1x,i6))
                    do 120 I=1,mul
                       if(FT5(I).le.t2time) then
                          gz_now=GZ5(i)+(t2time-FT5(i))*PZ5(i)/E5(i)
                       else
                          gz_now=GZ5(i)+(t2time-FT5(i))*vzp(i)
                       endif
                       If(dabs(gz_now).lt.t2time) then
                       etap=0.5d0*dlog((t2time+gz_now)/(t2time-gz_now))
                       else
                          etap=1000000.d0*sign(1.d0,gz_now)
                       endif
                       if(FT5(I).le.t2time) then
                          gx_now=GX5(i)+(t2time-FT5(i))*PX5(i)/E5(i)
                          gy_now=GY5(i)+(t2time-FT5(i))*PY5(i)/E5(i)
                       else
                          gx_now=GX5(i)+(t2time-FT5(i))*vxp(i)
                          gy_now=GY5(i)+(t2time-FT5(i))*vyp(i)
                       endif
                       write(93,185) ITYP5(i),PX5(i),PY5(i),PZ5(i),
     1                      XMASS5(i),gx_now,gy_now,ft5(i),etap
                       if(ifanim.eq.1.and.FT5(I).le.t2time) then
                          write(11,180) ITYP5(i),GX5(i),GY5(i),GZ5(i),
     1                         FT5(i),PX5(i),PY5(i),PZ5(i),E5(i)
                       endif
 185           format(i3,3(1x,f8.3),1x,f8.4,1x,2(f8.3,1x),f11.4,1x,f8.3)
 120                continue
                    itanim=nt
                 endif
 110          continue
           endif
 160       format(i10,3(2x,f9.5))
 180       format(i6,8(1x,f7.2))
        else if(idd.eq.3) then        
           do 1010 ianp=1,30
              if(iscatt(ianp).eq.0) tscatt(ianp+1)=tscatt(ianp)
 1010      continue
           do 350 I=1,mul
              delta=1d-8
              if((E5(i)-dabs(PZ5(i))+delta).le.0) then
                 write(6,*) 'ypartn error',E5(i)-dabs(PZ5(i))
                 ypartn=1000000.d0*sign(1.d0,PZ5(i))
              else
                 ypartn=0.5d0*dlog((E5(i)+PZ5(i)+delta)
     1                /(E5(i)-PZ5(i)+delta))
              endif
              pt2=PX5(I)**2+PY5(I)**2
              iloop=1
              if(dabs(ypartn).le.1d0) then
                 iloop=2
                 if(dabs(ypartn).le.0.5d0) then
                    iloop=3
                 endif
              endif
              do 325 ianp=1,30
                 if(iscatt(ianp).ne.0) then
                    if(FT5(I).lt.tscatt(ianp+1)
     1 .and.FT5(I).ge.tscatt(ianp)) then
                       do 1011 iy=1,iloop
                          if(pt2.gt.0d0) then
                             v2prtn=(PX5(I)**2-PY5(I)**2)/pt2
                             v2pfrz(ianp,iy)=v2pfrz(ianp,iy)+v2prtn
                     v2p2fz(ianp,iy)=v2p2fz(ianp,iy)+v2prtn**2
                          endif
                          xperp2=GX5(I)**2+GY5(I)**2
                          if(xperp2.gt.0d0) s2pfrz(ianp,iy)=
     1 s2pfrz(ianp,iy)+(GX5(I)**2-GY5(I)**2)/xperp2
        etpfrz(ianp,iy)=etpfrz(ianp,iy)+dsqrt(pt2+XMASS5(I)**2)
                          xnpfrz(ianp,iy)=xnpfrz(ianp,iy)+1d0
            if(ITYP5(I).eq.1.or.ITYP5(I).eq.2)then
              etlfrz(ianp,iy)=etlfrz(ianp,iy)+dsqrt(pt2+XMASS5(I)**2)
              xnlfrz(ianp,iy)=xnlfrz(ianp,iy)+1d0
            elseif(ITYP5(I).eq.3)then
              etsfrz(ianp,iy)=etsfrz(ianp,iy)+dsqrt(pt2+XMASS5(I)**2)
              xnsfrz(ianp,iy)=xnsfrz(ianp,iy)+1d0
            endif
 1011    continue
                       goto 350
                    endif
                 endif
 325          continue
 350       continue
        else if(idd.eq.2) then
           do 1012 i=1,3
              write(nfile(i),*) '   tsp,   v2p,     v2p2, '//
     1 '   s2p,  etp,   xmult,    nevt,  xnctot'
              write ((nfile(i)+1),*) '  tsp,   v2pf,   v2pf2, '//
     1 '   s2pf, etpf,  xnform,  xcrate'
              write ((nfile(i)+2),*) '  tsp,   v2pa,   v2pa2, '//
     1 '   s2pa, etpa,  xmult_ap,  xnform,   nevt'
              write ((nfile(i)+3),*) '  tsph,  v2ph,   v2ph2, '//
     1 '   s2ph, etph,  xmult_(ap/2+h),xmult_ap/2,nevt'
           write(nlfile(i),*) '   tsp,    v2,     s2,    etp,    xmul'
           write(nsfile(i),*) '   tsp,    v2,     s2,    etp,    xmul'
           write(nlfile(i)+1,*) '   tsp,    v2,     s2,    etp,    xmul'
           write(nsfile(i)+1,*) '   tsp,    v2,     s2,    etp,    xmul'
 1012   continue
           do 150 ii=1, 30
              if(nevt(ii).eq.0) goto 150
              do 1014 iy=1,3
                 if(xnpart(ii,iy).gt.1d0) then
                    v2p(ii,iy)=v2p(ii,iy)/xnpart(ii,iy)
                    v2p2(ii,iy)=dsqrt((v2p2(ii,iy)/xnpart(ii,iy)
     1                    -v2p(ii,iy)**2)/(xnpart(ii,iy)-1))
                    s2p(ii,iy)=s2p(ii,iy)/xnpart(ii,iy)
                    xmult=dble(xnpart(ii,iy)/dble(nevt(ii)))
                    etp(ii,iy)=etp(ii,iy)/dble(nevt(ii))
                    etpl(ii,iy)=etpl(ii,iy)/dble(nevt(ii))
                    etps(ii,iy)=etps(ii,iy)/dble(nevt(ii))
                    xnctot=0d0
                    do 1013 inum=1,ii
                       xnctot=xnctot+xncoll(inum)
 1013               continue
                    if(nevt(1).ne.0) xnctot=xnctot/nevt(1)
                    write (nfile(iy),200) tsp(ii),v2p(ii,iy),
     1      v2p2(ii,iy),s2p(ii,iy),etp(ii,iy),xmult,nevt(ii),xnctot
                 endif
                 if(nevt(ii).ne.0) 
     1                xcrate=xncoll(ii)/(tsp(ii+1)-tsp(ii))/nevt(ii)
                 if(xnpf(ii,iy).gt.1d0) then
                    v2pf(ii,iy)=v2pf(ii,iy)/xnpf(ii,iy)
                    v2pf2(ii,iy)=dsqrt((v2pf2(ii,iy)/xnpf(ii,iy)
     1                    -v2pf(ii,iy)**2)/(xnpf(ii,iy)-1))
                    s2pf(ii,iy)=s2pf(ii,iy)/xnpf(ii,iy)
                    xnform=dble(xnpf(ii,iy)/dble(nevt(ii)))
                    etpf(ii,iy)=etpf(ii,iy)/dble(nevt(ii))
                    etplf(ii,iy)=etplf(ii,iy)/dble(nevt(ii))
                    etpsf(ii,iy)=etpsf(ii,iy)/dble(nevt(ii))
                    write (nfile(iy)+1, 210) tsp(ii),v2pf(ii,iy),
     1      v2pf2(ii,iy),s2pf(ii,iy),etpf(ii,iy),xnform,xcrate
                 endif
                 if(xnpl(ii,iy).gt.1d0) then
                    v2pl(ii,iy)=v2pl(ii,iy)/xnpl(ii,iy)
                    s2pl(ii,iy)=s2pl(ii,iy)/xnpl(ii,iy)
                    xmult=dble(xnpl(ii,iy)/dble(nevt(ii)))
                    etpl(ii,iy)=etpl(ii,iy)/dble(nevt(ii))
                    write (nlfile(iy),201) tsp(ii),v2pl(ii,iy),
     1        s2pl(ii,iy),etpl(ii,iy),xmult
                 endif
                 if(xnps(ii,iy).gt.1d0) then
                    v2ps(ii,iy)=v2ps(ii,iy)/xnps(ii,iy)
                    s2ps(ii,iy)=s2ps(ii,iy)/xnps(ii,iy)
                    xmult=dble(xnps(ii,iy)/dble(nevt(ii)))
                    etps(ii,iy)=etps(ii,iy)/dble(nevt(ii))
                    write (nsfile(iy),201) tsp(ii),v2ps(ii,iy),
     1        s2ps(ii,iy),etps(ii,iy),xmult
                 endif
                 if(xnplf(ii,iy).gt.1d0) then
                    v2plf(ii,iy)=v2plf(ii,iy)/xnplf(ii,iy)
                    s2plf(ii,iy)=s2plf(ii,iy)/xnplf(ii,iy)
                    xmult=dble(xnplf(ii,iy)/dble(nevt(ii)))
                    etplf(ii,iy)=etplf(ii,iy)/dble(nevt(ii))
                    write (nlfile(iy)+1,201) tsp(ii),v2plf(ii,iy),
     1        s2plf(ii,iy),etplf(ii,iy),xmult
                 endif
                 if(xnpsf(ii,iy).gt.1d0) then
                    v2psf(ii,iy)=v2psf(ii,iy)/xnpsf(ii,iy)
                    s2psf(ii,iy)=s2psf(ii,iy)/xnpsf(ii,iy)
                    xmult=dble(xnpsf(ii,iy)/dble(nevt(ii)))
                    etpsf(ii,iy)=etpsf(ii,iy)/dble(nevt(ii))
                    write (nsfile(iy)+1,201) tsp(ii),v2psf(ii,iy),
     1        s2psf(ii,iy),etpsf(ii,iy),xmult
                 endif
 1014         continue
 150           continue
               scalei=0d0
               scalef=0d0               
               if(nevt(1).ne.0) SCALEi = 1d0 / dble(nevt(1)) / BMT
               if(nevt(30).ne.0) SCALEf = 1d0 / dble(nevt(30)) / BMT
         do 1016 iy=2,3
           yra = 1d0
           if(iy .eq. 2)yra = 2d0
         do 1015 i=1,50
           WRITE(nmfile(iy),251) BMT*dble(I - 0.5), 
     &     SCALEi*DMYil(I,iy)/yra, SCALEf*DMYfl(I,iy)/yra,
     &     SCALEi*DMYis(I,iy)/yra, SCALEf*DMYfs(I,iy)/yra
 1015   continue
 1016 continue
           if(nevt(30).ge.1) then
              do 1017 iy=1,3
                 v2pavg(iy)=v2psum(iy)/nevt(30)
                 v2var0=v2p2sm(iy)/nevt(30)-v2pavg(iy)**2
                 if(v2var0.gt.0d0) varv2p(iy)=dsqrt(v2var0)
 1017 continue
              write(49, 240) 'EBE v2p,v2p(y2),v2p(y1): avg=', v2pavg
              write(49, 240) 'EBE v2p,v2p(y2),v2p(y1): var=', varv2p
           endif
           if(ifanim.eq.1) then
              do 1018 I=1,mul
                 if(FT5(I).le.t2time) then
                    write(15,140) ITYP5(i),GX5(i),GY5(i),GZ5(i),FT5(i)
                 endif
 1018         continue
 140          format(i10,4(2x,f7.2))
              write(10,*) -10.
              write(10,*) 0
              write(11,*) -10.
              write(11,*) 0
              close(10)
              close(11)
              close(15)
           endif
           do 450 ianp=1,30
              do 400 iy=1,3
                 v2pact=0d0
                 v2p2ac=0d0
                 s2pact=0d0
                 etpact=0d0
                 xnacti=0d0
                 if(xnpf(ianp,iy).gt.1d0) then
                    v2pact=v2pf(ianp,iy)*xnpf(ianp,iy)
                    v2p2ac=(v2pf2(ianp,iy)**2*(xnpf(ianp,iy)-1)
     1 +v2pf(ianp,iy)**2)*xnpf(ianp,iy)
                    s2pact=s2pf(ianp,iy)*xnpf(ianp,iy)
                    etpact=etpf(ianp,iy)*dble(nevt(ianp))
                    xnpact=xnpf(ianp,iy)
                    do 1019 kanp=1,ianp
                       v2pact=v2pact-v2pfrz(kanp,iy)
                       v2p2ac=v2p2ac-v2p2fz(kanp,iy)
                       s2pact=s2pact-s2pfrz(kanp,iy)
                       etpact=etpact-etpfrz(kanp,iy)
                       xnpact=xnpact-xnpfrz(kanp,iy)
 1019               continue
                    v2ph=v2pact
                    v2ph2=v2p2ac
                    s2ph=s2pact
                    etph=etpact
                    xnp2=xnpact/2d0
                    if(xnpact.gt.1d0.and.nevt(ianp).ne.0) then
                       v2pact=v2pact/xnpact
                       v2p2ac=dsqrt((v2p2ac/xnpact
     1                    -v2pact**2)/(xnpact-1))
                       s2pact=s2pact/xnpact
                       xnacti=dble(xnpact/dble(nevt(ianp)))
                       etpact=etpact/dble(nevt(ianp))
                       write (nfile(iy)+2, 250) tsp(ianp),v2pact,
     1 v2p2ac,s2pact,etpact,xnacti,
     2 xnpf(ianp,iy)/dble(nevt(ianp)),nevt(ianp)
                    endif
                 endif
                 shadr=dble(nevt(ianp))/dble(nevent)
                 ianh=ianp
                 v2ph=v2ph+v2h(ianh,iy)*xnhadr(ianh,iy)*shadr
                 v2ph2=v2ph2+(v2h2(ianh,iy)**2*(xnhadr(ianh,iy)-1)
     1 +v2h(ianh,iy)**2)*xnhadr(ianh,iy)*shadr
                 s2ph=s2ph+s2h(ianh,iy)*xnhadr(ianh,iy)*shadr
                 etph=etph+eth(ianh,iy)*dble(nevent)*shadr
                 xnph=xnpact+xnhadr(ianh,iy)*shadr
                 xnp2h=xnp2+xnhadr(ianh,iy)*shadr
                 if(xnph.gt.1d0.and.nevt(ianp).ne.0) then
                    v2ph=v2ph/xnph
                    v2ph2=dsqrt((v2ph2/xnph-v2ph**2)/(xnph-1))
                    s2ph=s2ph/xnph
                    etph=etph/dble(nevt(ianp))
                    xnp2=xnp2/dble(nevt(ianp))
                    xnp2h=xnp2h/dble(nevent)
                    if(tsp(ianp).le.dble(ntmax*dt)) 
     1                    write (nfile(iy)+3, 250) tsp(ianp),v2ph,
     2 v2ph2,s2ph,etph,xnp2h,xnp2,nevt(ianp)
                 endif
 400              continue
 450       continue
           do 550 ianp=1,30
              do 500 iy=1,3
                 v2pact=0d0
                 v2p2ac=0d0
                 s2pact=0d0
                 etpact=0d0
                 xnacti=0d0
                    v2pact=v2pf(ianp,iy)*xnpf(ianp,iy)
                    v2p2ac=(v2pf2(ianp,iy)**2*(xnpf(ianp,iy)-1)
     1 +v2pf(ianp,iy)**2)*xnpf(ianp,iy)
                    s2pact=s2pf(ianp,iy)*xnpf(ianp,iy)
                    etpact=etpf(ianp,iy)*dble(nevt(ianp))
                    xnpact=xnpf(ianp,iy)
 500              continue
 550           continue
           close (620)
           close (630)
           do 1021 nf=1,3
              do 1020 ifile=0,3
                 close(nfile(nf)+ifile)
 1020        continue
 1021     continue
           do 1022 nf=1,3
              close(740+nf)
 1022      continue
        endif
 200        format(2x,f5.2,3(2x,f7.4),2(2x,f9.2),i6,2x,f9.2)
 210        format(2x,f5.2,3(2x,f7.4),3(2x,f9.2))
 240        format(a30,3(2x,f9.5))
 250        format(2x,f5.2,3(2x,f7.4),3(2x,f9.2),i6)
 201        format(2x,f5.2,4(2x,f9.2))
 251        format(5e15.5)
        return
        end
