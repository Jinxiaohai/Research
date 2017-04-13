        subroutine flowh0(NEVNT,idd)
        dimension tsh(31)
        DOUBLE PRECISION  v2h,xnhadr,eth,v2h2,s2h
        DOUBLE PRECISION  v2hp,xnhadp,v2hsum,v2h2sm,
     1 v2havg(3),varv2h(3)
        COMMON /hflow/ v2h(30,3),xnhadr(30,3),eth(30,3),
     1 v2h2(30,3),s2h(30,3)
        common/ebe/v2hp(3),xnhadp(3),v2hsum(3),v2h2sm(3)
        common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
        COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &   IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
        common /lastt/itimeh,bimp
        SAVE   
        if(idd.eq.0) then
           itimeh=0
           do 1001 ii = 1, 31
              tsh(ii)=float(ii-1)
 1001      continue
           do 1003 ii=1,30
              do 1002 iy=1,3
                 v2h(ii,iy)=0d0
                 xnhadr(ii,iy)=0d0
                 eth(ii,iy)=0d0
                 v2h2(ii,iy)=0d0
                 s2h(ii,iy)=0d0
 1002         continue
 1003      continue
           do 1004 iy=1,3
              v2hp(iy)=0d0
              xnhadp(iy)=0d0
              v2hsum(iy)=0d0
              v2h2sm(iy)=0d0
            if(iy.eq.1) then
               nunit=59
            elseif(iy.eq.2) then
               nunit=68
            else
               nunit=69
            endif
              write(nunit,*) '   tsh,   v2h,     v2h2,     s2h, '//
     1 ' eth,   xmulth'
 1004      continue
        else if(idd.eq.2) then
           do 100 ii=1, 30
              do 1005 iy=1,3
                 if(xnhadr(ii,iy).eq.0) then
                    xmulth=0.
                 elseif(xnhadr(ii,iy).gt.1) then
                    v2h(ii,iy)=v2h(ii,iy)/xnhadr(ii,iy)
                    eth(ii,iy)=eth(ii,iy)/dble(NEVNT)
                    v2h2(ii,iy)=dsqrt((v2h2(ii,iy)/xnhadr(ii,iy)
     1                    -v2h(ii,iy)**2)/(xnhadr(ii,iy)-1))
                    s2h(ii,iy)=s2h(ii,iy)/xnhadr(ii,iy)
                    xmulth=sngl(xnhadr(ii,iy)/NEVNT)
                 endif
             if(iy.eq.1) then
                nunit=59
             elseif(iy.eq.2) then
                nunit=68
             else
                nunit=69
             endif
                 if(tsh(ii).le.(ntmax*dt)) 
     1                    write (nunit,200) tsh(ii),v2h(ii,iy),
     2      v2h2(ii,iy),s2h(ii,iy),eth(ii,iy),xmulth
 1005         continue
 100           continue
           do 1006 iy=1,3
              v2havg(iy)=v2hsum(iy)/dble(NEVNT)
      varv2h(iy)=dsqrt(v2h2sm(iy)/dble(NEVNT)-v2havg(iy)**2)
 1006 continue
           write(88, 240) 'EBE v2h,v2h(y2),v2h(y1): avg=', v2havg
           write(88, 240) 'EBE v2h,v2h(y2),v2h(y1): var=', varv2h
        endif
 200        format(2x,f5.2,3(2x,f7.4),2(2x,f9.2))
 240        format(a30,3(2x,f9.5))
        return
        end
