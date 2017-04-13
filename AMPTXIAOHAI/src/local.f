        subroutine local(t)
        implicit double precision  (a-h, o-z)
        PARAMETER (MAXPTN=400001)
        PARAMETER (r0=1d0)
        COMMON /para1/ mul
        COMMON /prec2/GX5(MAXPTN),GY5(MAXPTN),GZ5(MAXPTN),FT5(MAXPTN),
     &       PX5(MAXPTN), PY5(MAXPTN), PZ5(MAXPTN), E5(MAXPTN),
     &       XMASS5(MAXPTN), ITYP5(MAXPTN)
        common /frzprc/ 
     &       gxfrz(MAXPTN), gyfrz(MAXPTN), gzfrz(MAXPTN), ftfrz(MAXPTN),
     &       pxfrz(MAXPTN), pyfrz(MAXPTN), pzfrz(MAXPTN), efrz(MAXPTN),
     &       xmfrz(MAXPTN), 
     &       tfrz(302), ifrz(MAXPTN), idfrz(MAXPTN), itlast
        common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
        common /prec5/ eta(MAXPTN), rap(MAXPTN), tau(MAXPTN)
        common /coal/dpcoal,drcoal,ecritl
        SAVE   
      do 1001 it=1,301
         if(t.ge.tfrz(it).and.t.lt.tfrz(it+1)) then
            if(it.eq.itlast) then
               return
            else
               itlast=it
               goto 50
            endif
         endif
 1001 continue
      write(1,*) 'local time out of range in LOCAL, stop',t,it
      stop
 50   continue
      do 200 ip=1,mul
         if(ifrz(ip).eq.1) goto 200
         if(it.eq.301) then
            etcrit=1d6
            goto 150
         else
            etcrit=(ecritl*2d0/3d0)
         endif
         if(t.lt.FT5(ip)) goto 200
         rap0=rap(ip)
         eta0=eta(ip)
         x0=GX5(ip)+vx(ip)*(t-FT5(ip))
         y0=GY5(ip)+vy(ip)*(t-FT5(ip))
         detdy=0d0
         do 100 itest=1,mul
            if(itest.eq.ip.or.t.lt.FT5(itest)) goto 100
            ettest=eta(itest)
            xtest=GX5(itest)+vx(itest)*(t-FT5(itest))
            ytest=GY5(itest)+vy(itest)*(t-FT5(itest))
            drt=sqrt((xtest-x0)**2+(ytest-y0)**2)
            if(dabs(ettest-eta0).le.1d0.and.drt.le.r0) 
     1           detdy=detdy+dsqrt(PX5(itest)**2+PY5(itest)**2
     2           +XMASS5(itest)**2)*0.5d0
 100     continue
         detdy=detdy*(dcosh(eta0)**2)/(t*3.1416d0*r0**2*dcosh(rap0))
 150     if(detdy.le.etcrit) then
            ifrz(ip)=1
            idfrz(ip)=ITYP5(ip)
            pxfrz(ip)=PX5(ip)
            pyfrz(ip)=PY5(ip)
            pzfrz(ip)=PZ5(ip)
            efrz(ip)=E5(ip)
            xmfrz(ip)=XMASS5(ip)
            if(t.gt.FT5(ip)) then
               gxfrz(ip)=x0
               gyfrz(ip)=y0
               gzfrz(ip)=GZ5(ip)+vz(ip)*(t-FT5(ip))
               ftfrz(ip)=t
            else
               gxfrz(ip)=GX5(ip)
               gyfrz(ip)=GY5(ip)
               gzfrz(ip)=GZ5(ip)
               ftfrz(ip)=FT5(ip)
            endif
         endif
 200  continue
        return
        end
