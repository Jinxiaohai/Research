      subroutine hbtout(nnew,nt,ntmax)
      PARAMETER  (MAXSTR=150001,MAXR=1)
      PARAMETER  (oneminus=0.99999,oneplus=1.00001)
      dimension lastkp(MAXSTR), newkp(MAXSTR),xnew(3)
      common /para7/ ioscar,nsmbbbar,nsmmeson
      COMMON/hbt/lblast(MAXSTR),xlast(4,MAXSTR),plast(4,MAXSTR),nlast
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      COMMON   /AA/  R(3,MAXSTR)
      COMMON   /BB/  P(3,MAXSTR)
      COMMON   /CC/  E(MAXSTR)
      COMMON   /EE/  ID(MAXSTR),LB(MAXSTR)
      common /lastt/itimeh,bimp
      COMMON/tdecay/tfdcy(MAXSTR),tfdpi(MAXSTR,MAXR),tft(MAXSTR)
      COMMON /AREVT/ IAEVT, IARUN, MISS
      common/snn/efrm,npart1,npart2,epsiPz,epsiPt,PZPROJ,PZTARG
      COMMON/HJGLBR/NELT,NINTHJ,NELP,NINP
      COMMON/FTMAX/ftsv(MAXSTR),ftsvt(MAXSTR, MAXR)
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      EXTERNAL IARFLV, INVFLV
      common /para8/ idpert,npertd,idxsec
      common /phiHJ/iphirp,phiRP
      SAVE   
      do 1001 i=1,max0(nlast,nnew)
         lastkp(i)=0
 1001 continue
      do 1002 i=1,nnew
         newkp(i)=0
 1002 continue
      do 100 ip=1,nnew
         do 1004 iplast=1,nlast
            if(p(1,ip).eq.plast(1,iplast).and.
     1           p(2,ip).eq.plast(2,iplast).and.
     2           p(3,ip).eq.plast(3,iplast).and.
     3           e(ip).eq.plast(4,iplast).and.
     4           lb(ip).eq.lblast(iplast).and.
     5      dpertp(ip).eq.dplast(iplast).and.lastkp(iplast).eq.0) then
               deltat=nt*dt-xlast(4,iplast)
               ene=sqrt(plast(1,iplast)**2+plast(2,iplast)**2
     1              +plast(3,iplast)**2+plast(4,iplast)**2)
               do 1003 ii=1,3
                  xnew(ii)=xlast(ii,iplast)+plast(ii,iplast)/ene*deltat
 1003          continue
                  dr=sqrt((r(1,ip)-xnew(1))**2+(r(2,ip)-xnew(2))**2
     1              +(r(3,ip)-xnew(3))**2)
               if(dr.le.0.01) then
                  lastkp(iplast)=1
                  newkp(ip)=1
                  if(nt.eq.ntmax.and.ftsv(ip).gt.((ntmax-1)*dt)) 
     1                 xlast(4,iplast)=ftsv(ip)
                  goto 100
               endif
            endif
 1004    continue
 100  continue
      do 150 ip=1,nnew
         if(newkp(ip).eq.0) then
            do 1005 iplast=1,nnew
               if(lastkp(iplast).eq.0) then
                  xlast(1,iplast)=r(1,ip)
                  xlast(2,iplast)=r(2,ip)
                  xlast(3,iplast)=r(3,ip)
                  xlast(4,iplast)=nt*dt
                  if(nt.eq.ntmax) then
                     if(tfdcy(ip).gt.(ntmax*dt+0.001)) then
                        xlast(4,iplast)=tfdcy(ip)
                     elseif(ftsv(ip).gt.((ntmax-1)*dt)) then
                        xlast(4,iplast)=ftsv(ip)
                     endif
                  endif
                  plast(1,iplast)=p(1,ip)
                  plast(2,iplast)=p(2,ip)
                  plast(3,iplast)=p(3,ip)
                  plast(4,iplast)=e(ip)
                  lblast(iplast)=lb(ip)
                  lastkp(iplast)=1
                  dplast(iplast)=dpertp(ip)
                  goto 150
               endif
 1005       continue
         endif
 150  continue
      if(nnew.lt.nlast) then
         do 170 iplast=1,nlast
            if(lastkp(iplast).eq.0) then
               do 1006 ip2=iplast+1,nlast
                  if(lastkp(ip2).eq.1) then
                     xlast(1,iplast)=xlast(1,ip2)
                     xlast(2,iplast)=xlast(2,ip2)
                     xlast(3,iplast)=xlast(3,ip2)
                     xlast(4,iplast)=xlast(4,ip2)
                     plast(1,iplast)=plast(1,ip2)
                     plast(2,iplast)=plast(2,ip2)
                     plast(3,iplast)=plast(3,ip2)
                     plast(4,iplast)=plast(4,ip2)
                     lblast(iplast)=lblast(ip2)
                     lastkp(iplast)=1
                     dplast(iplast)=dplast(ip2)
                     goto 170
                  endif
 1006          continue
            endif
 170     continue
      endif
      nlast=nnew
      if(nt.eq.ntmax) then
         ndpert=0
         do ip=1,nlast
            if(dplast(ip).gt.oneminus.and.dplast(ip).lt.oneplus) then
            else
               ndpert=ndpert+1
            endif
         enddo
         write(16,191) IAEVT,IARUN,nlast-ndpert,bimp,npart1,npart2,
     1 NELP,NINP,NELT,NINTHJ,phiRP
         if(idpert.eq.1.or.idpert.eq.2)
     1        write(90,190) IAEVT,IARUN,ndpert,bimp,npart1,npart2,
     2        NELP,NINP,NELT,NINTHJ
         do 1007 ip=1,nlast
            if(abs(plast(1,ip)).le.epsiPt.and.abs(plast(2,ip)).le.epsiPt
     1           .and.(plast(3,ip).gt.amax1(0.,PZPROJ-epsiPz)
     2                .or.plast(3,ip).lt.(-PZTARG+epsiPz))
     3           .and.(lblast(ip).eq.1.or.lblast(ip).eq.2)) then
               if(dplast(ip).gt.oneminus.and.dplast(ip).lt.oneplus) then
                  write(16,200) INVFLV(lblast(ip)), plast(1,ip),
     1                 plast(2,ip),plast(3,ip),plast(4,ip),
     2                 xlast(1,ip),xlast(2,ip),xlast(3,ip),
     3                 xlast(4,ip)
               else
                  if(idpert.eq.1.or.idpert.eq.2) then
                     write(90,250) INVFLV(lblast(ip)), plast(1,ip),
     1                 plast(2,ip),plast(3,ip),
     2                 xlast(1,ip),xlast(2,ip),xlast(3,ip),
     3                 xlast(4,ip)
                  else
                     write(99,*) 'Unexpected perturbative particles'
                  endif
               endif
            elseif(amax1(abs(xlast(1,ip)),abs(xlast(2,ip)),
     1              abs(xlast(3,ip)),abs(xlast(4,ip))).lt.9999) then
               if(dplast(ip).gt.oneminus.and.dplast(ip).lt.oneplus) then
            write(16,200) INVFLV(lblast(ip)), plast(1,ip),
     1           plast(2,ip),plast(3,ip),plast(4,ip),
     2           xlast(1,ip),xlast(2,ip),xlast(3,ip),xlast(4,ip)
               else
                  if(idpert.eq.1.or.idpert.eq.2) then
            write(90,250) INVFLV(lblast(ip)),plast(1,ip),
     1           plast(2,ip),plast(3,ip),
     2           xlast(1,ip),xlast(2,ip),xlast(3,ip),xlast(4,ip),
     3           dplast(ip)
                  else
                     write(99,*) 'Unexpected perturbative particles'
                  endif
               endif
            else
               if(dplast(ip).gt.oneminus.and.dplast(ip).lt.oneplus) then
            write(16,201) INVFLV(lblast(ip)), plast(1,ip),
     1           plast(2,ip),plast(3,ip),plast(4,ip),
     2           xlast(1,ip),xlast(2,ip),xlast(3,ip),xlast(4,ip)
               else
                  if(idpert.eq.1.or.idpert.eq.2) then
                     write(90,251) INVFLV(lblast(ip)), plast(1,ip),
     1           plast(2,ip),plast(3,ip),
     2           xlast(1,ip),xlast(2,ip),xlast(3,ip),xlast(4,ip),
     3           dplast(ip)
                  else
                     write(99,*) 'Unexpected perturbative particles'
                  endif
               endif
            endif
 1007    continue
         if(ioscar.eq.1) call hoscar
      endif
 190  format(3(i7),f10.4,5x,6(i4))
 191  format(3(i7),f10.4,5x,6(i4),5x,f7.4)
 200  format(I6,2(1x,f8.3),1x,f11.4,1x,f6.3,4(1x,f8.2))
 201  format(I6,2(1x,f8.3),1x,f11.4,1x,f6.3,4(1x,e8.2))
 250  format(I5,2(1x,f8.3),1x,f10.3,2(1x,f7.1),1x,f8.2,1x,f7.2,1x,e10.4)
 251  format(I5,2(1x,f8.3),1x,f10.3,4(1x,e8.2),1x,e10.4)
        return
        end
