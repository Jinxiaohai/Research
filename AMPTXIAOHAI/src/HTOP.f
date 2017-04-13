      SUBROUTINE HTOP
      PARAMETER (MAXSTR=150001)
      PARAMETER (MAXPTN=400001)
      PARAMETER (MAXIDL=4001)
      DOUBLE PRECISION  GX0, GY0, GZ0, FT0, PX0, PY0, PZ0, E0, XMASS0
      DOUBLE PRECISION  PXSGS,PYSGS,PZSGS,PESGS,PMSGS,
     1     GXSGS,GYSGS,GZSGS,FTSGS, ptwo, xmdq, ptwox, ptwoy, ptwoz
      dimension it(4)
      COMMON/HMAIN2/KATT(MAXSTR,4),PATT(MAXSTR,4)
      COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
      COMMON /PARA1/ MUL
      COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &     PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &     XMASS0(MAXPTN), ITYP0(MAXPTN)
      COMMON /ilist7/ LSTRG0(MAXPTN), LPART0(MAXPTN)
      COMMON /ARPRC/ ITYPAR(MAXSTR),
     &     GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &     PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &     XMAR(MAXSTR)
      common /decom/ptwo(2,5)
      COMMON/RNDF77/NSEED
      COMMON /NOPREC/ NNOZPC, ITYPN(MAXIDL),
     &     GXN(MAXIDL), GYN(MAXIDL), GZN(MAXIDL), FTN(MAXIDL),
     &     PXN(MAXIDL), PYN(MAXIDL), PZN(MAXIDL), EEN(MAXIDL),
     &     XMN(MAXIDL)
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      COMMON/SOFT/PXSGS(MAXSTR,3),PYSGS(MAXSTR,3),PZSGS(MAXSTR,3),
     &     PESGS(MAXSTR,3),PMSGS(MAXSTR,3),GXSGS(MAXSTR,3),
     &     GYSGS(MAXSTR,3),GZSGS(MAXSTR,3),FTSGS(MAXSTR,3),
     &     K1SGS(MAXSTR,3),K2SGS(MAXSTR,3),NJSGS(MAXSTR)
      common/anim/nevent,isoft,isflag,izpc
      DOUBLE PRECISION vxp0,vyp0,vzp0,xstrg0,ystrg0,xstrg,ystrg
      common /precpa/vxp0(MAXPTN),vyp0(MAXPTN),vzp0(MAXPTN),
     1     xstrg0(MAXPTN),ystrg0(MAXPTN),
     2     xstrg(MAXPTN),ystrg(MAXPTN),istrg0(MAXPTN),istrg(MAXPTN)
      common /para7/ ioscar,nsmbbbar,nsmmeson
      COMMON /AREVT/ IAEVT, IARUN, MISS
      common/snn/efrm,npart1,npart2,epsiPz,epsiPt,PZPROJ,PZTARG
      SAVE   
        npar=0
        nnozpc=0
        if((isoft.eq.4.or.isoft.eq.5).and.(ioscar.eq.2.or.ioscar.eq.3)) 
     1       then
           nsmbbbar=0
           nsmmeson=0
           do i=1,natt
              id=ITYPAR(i)
              idabs=iabs(id)
              i2=MOD(idabs/10,10)
              if(abs(PXAR(i)).le.epsiPt.and.abs(PYAR(i)).le.epsiPt
     1             .and.(PZAR(i).gt.amax1(0.,PZPROJ-epsiPz)
     2                .or.PZAR(i).lt.(-PZTARG+epsiPz))
     3             .and.(id.eq.2112.or.id.eq.2212)) then
              elseif(idabs.gt.1000.and.i2.ne.0) then
                 nsmbbbar=nsmbbbar+1
              elseif((idabs.gt.100.and.idabs.lt.1000)
     1                .or.idabs.gt.10000) then
                 nsmmeson=nsmmeson+1
              endif
           enddo
           if(ioscar.eq.2.or.ioscar.eq.3) then
              write(92,*) iaevt,miss,3*nsmbbbar+2*nsmmeson,
     1             nsmbbbar,nsmmeson,natt,natt-nsmbbbar-nsmmeson
           endif
        endif
        do 100 i=1,natt
           id=ITYPAR(i)
           idabs=iabs(id)
           i4=MOD(idabs/1000,10)
           i3=MOD(idabs/100,10)
           i2=MOD(idabs/10,10)
           i1=MOD(idabs,10)
           rnum=RANART(NSEED)
           ftime=0.197*PEAR(i)/(PXAR(i)**2+PYAR(i)**2+XMAR(i)**2)
           inozpc=0
           it(1)=0
           it(2)=0
           it(3)=0
           it(4)=0
              if(abs(PXAR(i)).le.epsiPt.and.abs(PYAR(i)).le.epsiPt
     1             .and.(PZAR(i).gt.amax1(0.,PZPROJ-epsiPz)
     2                .or.PZAR(i).lt.(-PZTARG+epsiPz))
     3             .and.(id.eq.2112.or.id.eq.2212)) then
              inozpc=1
           elseif(idabs.gt.1000.and.i2.ne.0) then
              if(((i4.eq.1.or.i4.eq.2).and.i4.eq.i3)
     1 .or.(i4.eq.3.and.i3.eq.3)) then
                 if(i1.eq.2) then
                    if(rnum.le.(1./2.)) then
                       it(1)=i4
                       it(2)=i3*1000+i2*100+1
                    elseif(rnum.le.(2./3.)) then
                       it(1)=i4
                       it(2)=i3*1000+i2*100+3
                    else
                       it(1)=i2
                       it(2)=i4*1000+i3*100+3
                    endif
                 elseif(i1.eq.4) then
                    if(rnum.le.(2./3.)) then
                       it(1)=i4
                       it(2)=i3*1000+i2*100+3
                    else
                       it(1)=i2
                       it(2)=i4*1000+i3*100+3
                    endif
                 endif
              elseif(i4.eq.1.or.i4.eq.2) then
                 if(i1.eq.2) then
                    if(rnum.le.(1./2.)) then
                       it(1)=i2
                       it(2)=i4*1000+i3*100+1
                    elseif(rnum.le.(2./3.)) then
                       it(1)=i2
                       it(2)=i4*1000+i3*100+3
                    else
                       it(1)=i4
                       it(2)=i3*1000+i2*100+3
                    endif
                 elseif(i1.eq.4) then
                    if(rnum.le.(2./3.)) then
                       it(1)=i2
                       it(2)=i4*1000+i3*100+3
                    else
                       it(1)=i4
                       it(2)=i3*1000+i2*100+3
                    endif
                 endif
              elseif(i4.ge.3) then
                 it(1)=i4
                 if(i3.lt.i2) then
                    it(2)=i2*1000+i3*100+1
                 else
                    it(2)=i3*1000+i2*100+3
                 endif
              endif
              if(id.lt.0) then
                 it(1)=-it(1)
                 it(2)=-it(2)
              endif
              if(isoft.eq.4.or.isoft.eq.5) then
                 it(3)=MOD(it(2)/1000,10)
                 it(4)=MOD(it(2)/100,10)
              endif
           elseif((idabs.gt.100.and.idabs.lt.1000)
     1 .or.idabs.gt.10000) then
              if(i3.eq.i2) then
                 if(i3.eq.1.or.i3.eq.2) then
                    if(rnum.le.0.5) then
                       it(1)=1
                       it(2)=-1
                    else
                       it(1)=2
                       it(2)=-2
                    endif
                 else
                    it(1)=i3
                    it(2)=-i3
                 endif
              else
                 if((isign(1,id)*(-1)**i3).eq.1) then
                    it(1)=i3
                    it(2)=-i2
                 else
                    it(1)=i2
                    it(2)=-i3
                 endif
              endif
           else
              inozpc=1
           endif
           if(inozpc.eq.1) then
              NJSGS(i)=0
              nnozpc=nnozpc+1
              itypn(nnozpc)=ITYPAR(i)
              pxn(nnozpc)=PXAR(i)
              pyn(nnozpc)=PYAR(i)
              pzn(nnozpc)=PZAR(i)
              een(nnozpc)=PEAR(i)
              xmn(nnozpc)=XMAR(i)
              gxn(nnozpc)=GXAR(i)
              gyn(nnozpc)=GYAR(i)
              gzn(nnozpc)=GZAR(i)
              ftn(nnozpc)=FTAR(i)
           else
              NJSGS(i)=2
              ptwo(1,5)=dble(ulmass(it(1)))
              ptwo(2,5)=dble(ulmass(it(2)))
              call decomp(dble(patt(i,1)),dble(patt(i,2)),
     1             dble(patt(i,3)),dble(XMAR(i)),i,it(1))
              ipamax=2
              if((isoft.eq.4.or.isoft.eq.5)
     1 .and.iabs(it(2)).gt.1000) ipamax=1
              do 1001 ipar=1,ipamax
                 npar=npar+1
                 ityp0(npar)=it(ipar)
                 px0(npar)=ptwo(ipar,1)
                 py0(npar)=ptwo(ipar,2)
                 pz0(npar)=ptwo(ipar,3)
                 e0(npar)=ptwo(ipar,4)
                 xmass0(npar)=ptwo(ipar,5)
                 gx0(npar)=dble(GXAR(i))
                 gy0(npar)=dble(GYAR(i))
                 gz0(npar)=dble(GZAR(i))
                 ft0(npar)=dble(ftime)
                 lstrg0(npar)=i
                 lpart0(npar)=ipar
                 vxp0(npar)=dble(patt(i,1)/patt(i,4))
                 vyp0(npar)=dble(patt(i,2)/patt(i,4))
                 vzp0(npar)=dble(patt(i,3)/patt(i,4))
                 xstrg(npar)=xstrg0(i)
                 ystrg(npar)=ystrg0(i)
                 istrg(npar)=istrg0(i)
 1001     continue
 200      format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f8.2))
 201      format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e8.2))
              if((isoft.eq.4.or.isoft.eq.5)
     1 .and.iabs(it(2)).gt.1000) then
                 NJSGS(i)=3
                 xmdq=ptwo(2,5)
                 ptwo(1,5)=dble(ulmass(it(3)))
                 ptwo(2,5)=dble(ulmass(it(4)))
                 ptwox=ptwo(2,1)
                 ptwoy=ptwo(2,2)
                 ptwoz=ptwo(2,3)
                 call decomp(ptwox,ptwoy,ptwoz,xmdq,i,it(1))
                 do 1002 ipar=1,2
                    npar=npar+1
                    ityp0(npar)=it(ipar+2)
                    px0(npar)=ptwo(ipar,1)
                    py0(npar)=ptwo(ipar,2)
                    pz0(npar)=ptwo(ipar,3)
                    e0(npar)=ptwo(ipar,4)
                    xmass0(npar)=ptwo(ipar,5)
                    gx0(npar)=dble(GXAR(i))
                    gy0(npar)=dble(GYAR(i))
                    gz0(npar)=dble(GZAR(i))
                    ft0(npar)=dble(ftime)
                    lstrg0(npar)=i
                    lpart0(npar)=ipar+1
                    vxp0(npar)=dble(patt(i,1)/patt(i,4))
                    vyp0(npar)=dble(patt(i,2)/patt(i,4))
                    vzp0(npar)=dble(patt(i,3)/patt(i,4))
                    xstrg(npar)=xstrg0(i)
                    ystrg(npar)=ystrg0(i)
                    istrg(npar)=istrg0(i)
 1002        continue
              endif
           endif
 100        continue
      MUL=NPAR
      if((isoft.eq.4.or.isoft.eq.5).and.(ioscar.eq.2.or.ioscar.eq.3)) 
     1     then
         if((natt-nsmbbbar-nsmmeson).ne.nnozpc) 
     1        write(92,*) 'Problem with the total # of initial particles
     2 (gamma,e,muon,...) not entering ZPC'
         if((3*nsmbbbar+2*nsmmeson).ne.npar) 
     1        write(92,*) 'Problem with the total # of initial partons
     2 after string melting'
      endif
      RETURN
      END
