        subroutine flowh(ct)
        PARAMETER (MAXSTR=150001, MAXR=1)
        dimension tsh(31)
        DOUBLE PRECISION  v2h,xnhadr,eth,v2h2,s2h
        DOUBLE PRECISION  v2hp,xnhadp,v2hsum,v2h2sm,v2hevt(3)
        DOUBLE PRECISION  pt2, v2hadr
        COMMON /hflow/ v2h(30,3),xnhadr(30,3),eth(30,3),
     1 v2h2(30,3),s2h(30,3)
        common/ebe/v2hp(3),xnhadp(3),v2hsum(3),v2h2sm(3)
        common /lastt/itimeh,bimp
        COMMON /RUN/ NUM
        COMMON  /AA/      R(3,MAXSTR)
        COMMON  /BB/      P(3,MAXSTR)
        COMMON  /CC/      E(MAXSTR)
        COMMON  /EE/      ID(MAXSTR),LB(MAXSTR)
        COMMON  /RR/      MASSR(0:MAXR)
        common/anim/nevent,isoft,isflag,izpc
        COMMON /AREVT/ IAEVT, IARUN, MISS
        SAVE   
        do 1001 ii = 1, 31
           tsh(ii)=float(ii-1)
 1001   continue
        do 1004 ianh = 1, 30
           if ((ct+0.0001).lt.tsh(ianh+1)
     1 .and.(ct+0.0001).ge.tsh(ianh)) then
              if(ianh.eq.itimeh) goto 101
              IA = 0
              DO 1002 J = 1, NUM
                 mult=MASSR(J)
                 IA = IA + MASSR(J - 1)
                 DO 100 IC = 1, mult
                    I = IA + IC
                    if(iabs(LB(I)-10000).lt.100) goto 100
                    px=p(1,i)
                    py=p(2,i)
                    pt2=dble(px)**2+dble(py)**2
                    ene=sqrt(e(i)**2+sngl(pt2)+p(3,i)**2)
                    RAP=0.5*alog((ene+p(3,i))/(ene-p(3,i)))
                    iloop=1
                    if(abs(rap).le.1) then
                       iloop=2
                       if(abs(rap).le.0.5) then
                          iloop=3
                       endif
                    endif
                    do 50 iy=1,iloop
                       if(pt2.gt.0d0) then
                          v2hadr=(dble(px)**2-dble(py)**2)/pt2
                          v2h(ianh,iy)=v2h(ianh,iy)+v2hadr
                          v2h2(ianh,iy)=v2h2(ianh,iy)+v2hadr**2
                          if(dabs(v2hadr).gt.1d0) 
     1 write(1,*) 'v2hadr>1',v2hadr,px,py
                       endif
                       xperp2=r(1,I)**2+r(2,I)**2
                       if(xperp2.gt.0.) 
     1 s2h(ianh,iy)=s2h(ianh,iy)+dble((r(1,I)**2-r(2,I)**2)/xperp2)
               eth(ianh,iy)=eth(ianh,iy)+dble(SQRT(e(i)**2+sngl(pt2)))
                       xnhadr(ianh,iy)=xnhadr(ianh,iy)+1d0
 50                    continue
 100                 continue
 1002         CONTINUE
              itimeh=ianh
              if(ianh.eq.30) then
                 do 1003 iy=1,3
                    nhadrn=IDINT(xnhadr(ianh,iy)-xnhadp(iy))
                    if(nhadrn.ne.0) then
               v2hevt(iy)=(v2h(ianh,iy)-v2hp(iy))/dble(nhadrn)
                       v2hsum(iy)=v2hsum(iy)+v2hevt(iy)
                       v2h2sm(iy)=v2h2sm(iy)+v2hevt(iy)**2
                       v2hp(iy)=v2h(ianh,iy)
                       xnhadp(iy)=xnhadr(ianh,iy)
                    endif
 1003            continue
                 write(88, 160) iaevt,v2hevt
              endif
              goto 101
           endif
 1004   continue
 160    format(i10,3(2x,f9.5))
 101    ifanim=0
        if(ifanim.eq.1) then
           IA = 0
           do 1005 J = 1, NUM
              mult=MASSR(J)
              IA = IA + MASSR(J - 1)
              write(10,*) ct
              write(10,*) mult
              DO 150 IC = 1, mult
                 I = IA + IC
                 if(amax1(abs(r(1,i)),abs(r(2,i)),
     1                abs(r(3,i))).lt.9999) then
                    write(10,210) LB(i),r(1,i),r(2,i),r(3,i),
     1                   p(1,i),p(2,i),p(3,i),e(i)
                 else
                    write(10,220) LB(i),r(1,i),r(2,i),r(3,i),
     1                   p(1,i),p(2,i),p(3,i),e(i)
                 endif
 150          continue
 1005      continue
           return
        endif
 210    format(i6,7(1x,f9.3))
 220    format(i6,3(1x,e9.3),4(1x,f9.3))
        return
        end
