      SUBROUTINE PTOH
      PARAMETER (MAXSTR=150001)
      DOUBLE PRECISION  gxp,gyp,gzp,ftp,pxp,pyp,pzp,pep,pmp
      DOUBLE PRECISION  gxp0,gyp0,gzp0,ft0fom,drlocl
      DOUBLE PRECISION  enenew, pxnew, pynew, pznew, beta2, gam
      DOUBLE PRECISION  ftavg0,gxavg0,gyavg0,gzavg0,bex,bey,bez
      DOUBLE PRECISION  PXSGS,PYSGS,PZSGS,PESGS,PMSGS,
     1     GXSGS,GYSGS,GZSGS,FTSGS
      DOUBLE PRECISION  xmdiag,px1,py1,pz1,e1,px2,py2,pz2,e2,
     1     px3,py3,pz3,e3,xmpair,etot
      DOUBLE PRECISION  p1,p2,p3
      common /loclco/gxp(3),gyp(3),gzp(3),ftp(3),
     1     pxp(3),pyp(3),pzp(3),pep(3),pmp(3)
      COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
      COMMON/HMAIN2/KATT(MAXSTR,4),PATT(MAXSTR,4)
      COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &     K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &     PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)
      COMMON /ARPRC/ ITYPAR(MAXSTR),
     &     GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &     PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &     XMAR(MAXSTR)
      COMMON/SOFT/PXSGS(MAXSTR,3),PYSGS(MAXSTR,3),PZSGS(MAXSTR,3),
     &     PESGS(MAXSTR,3),PMSGS(MAXSTR,3),GXSGS(MAXSTR,3),
     &     GYSGS(MAXSTR,3),GZSGS(MAXSTR,3),FTSGS(MAXSTR,3),
     &     K1SGS(MAXSTR,3),K2SGS(MAXSTR,3),NJSGS(MAXSTR)
      COMMON/RNDF77/NSEED
      common/anim/nevent,isoft,isflag,izpc
      common /prtn23/ gxp0(3),gyp0(3),gzp0(3),ft0fom
      common /nzpc/nattzp
      common /lor/ enenew, pxnew, pynew, pznew
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      common /lastt/itimeh,bimp
      COMMON/HJGLBR/NELT,NINTHJ,NELP,NINP
      COMMON /AREVT/ IAEVT, IARUN, MISS
      common /para7/ ioscar,nsmbbbar,nsmmeson
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      dimension xmdiag(MAXSTR),indx(MAXSTR),ndiag(MAXSTR)
      SAVE   
      call coales
      mstj24=MSTJ(24)
      MSTJ(24)=0
        nuudd=0
        npich=0
        nrhoch=0
      ppi0=1.
      prho0=0.
      DO 1001 ISG = 1, NSG
           if(NJSGS(ISG).ne.0) then
              NATT=NATT+1
              K1=K2SGS(ISG,1)
              k1abs=iabs(k1)
              PX1=PXSGS(ISG,1)
              PY1=PYSGS(ISG,1)
              PZ1=PZSGS(ISG,1)
              K2=K2SGS(ISG,2)
              k2abs=iabs(k2)
              PX2=PXSGS(ISG,2)
              PY2=PYSGS(ISG,2)
              PZ2=PZSGS(ISG,2)
              e1=PESGS(ISG,1)
              e2=PESGS(ISG,2)
              xmpair=dsqrt((e1+e2)**2-(px1+px2)**2-(py1+py2)**2
     1 -(pz1+pz2)**2)
              ibs=2
              imspin=0
              if(k1.eq.-k2.and.iabs(k1).le.2.
     1           and.NJSGS(ISG).eq.2) then
               nuudd=nuudd+1
               xmdiag(nuudd)=xmpair
               ndiag(nuudd)=natt
            endif
              K3=0
              if((isoft.eq.4.or.isoft.eq.5).and.NJSGS(ISG).eq.3) then
               K3=K2SGS(ISG,3)
               k3abs=iabs(k3)
               PX3=PXSGS(ISG,3)
               PY3=PYSGS(ISG,3)
               PZ3=PZSGS(ISG,3)
               e3=PESGS(ISG,3)
               xmpair=dsqrt((e1+e2+e3)**2-(px1+px2+px3)**2
     1              -(py1+py2+py3)**2-(pz1+pz2+pz3)**2)
              endif
              if(isoft.eq.3.and.
     1           (k1abs.gt.1000.or.k2abs.gt.1000)) then
               if(k1abs.gt.1000) then
                  kdq=k1abs
                  kk=k2abs
               else
                  kdq=k2abs
                  kk=k1abs
               endif
               ki=MOD(kdq/1000,10)
               kj=MOD(kdq/100,10)
               if(MOD(kdq,10).eq.1) then
                  idqspn=0
               else
                  idqspn=1
               endif
               if(kk.gt.ki) then
                  ktemp=kk
                  kk=kj
                  kj=ki
                  ki=ktemp
               elseif(kk.gt.kj) then
                  ktemp=kk
                  kk=kj
                  kj=ktemp
               endif
               if(ki.ne.kj.and.ki.ne.kk.and.kj.ne.kk) then
                  if(idqspn.eq.0) then
                     kf=1000*ki+100*kk+10*kj+ibs
                  else
                     kf=1000*ki+100*kj+10*kk+ibs
                  endif
               elseif(ki.eq.kj.and.ki.eq.kk) then
                  kf=1000*ki+100*kj+10*kk+4
               else
                  kf=1000*ki+100*kj+10*kk+ibs
               endif
               if(kf.eq.2112.or.kf.eq.2212) then
                  if(abs(sngl(xmpair)-ULMASS(kf)).gt.
     1                 abs(sngl(xmpair)-ULMASS(kf+2))) kf=kf+2
               endif
               if(k1.lt.0) kf=-kf
              elseif((isoft.eq.4.or.isoft.eq.5).and.NJSGS(ISG).eq.3) 
     1              then
               if(k1abs.gt.k2abs) then
                  ki=k1abs
                  kk=k2abs
               else
                  ki=k2abs
                  kk=k1abs
               endif
               if(k3abs.gt.ki) then
                  kj=ki
                  ki=k3abs
               elseif(k3abs.lt.kk) then
                  kj=kk
                  kk=k3abs
               else
                  kj=k3abs
               endif
               if(ki.eq.kj.and.ki.eq.kk) then
                  ibs=4
                  kf=1000*ki+100*kj+10*kk+ibs
               elseif(ki.ne.kj.and.ki.ne.kk.and.kj.ne.kk) then
                  ibs=2
                  kf1=1000*ki+100*kj+10*kk+ibs
                  kf2=1000*ki+100*kk+10*kj+ibs
                  kf=kf1
                  if(abs(sngl(xmpair)-ULMASS(kf1)).gt.
     1                 abs(sngl(xmpair)-ULMASS(kf2))) kf=kf2
               else
                  ibs=2
                  kf=1000*ki+100*kj+10*kk+ibs
                  if(kf.eq.2112.or.kf.eq.2212) then
                     if(abs(sngl(xmpair)-ULMASS(kf)).gt.
     1                    abs(sngl(xmpair)-ULMASS(kf+2))) kf=kf+2
                  endif
               endif
               if(k1.lt.0) kf=-kf
              else
               if(k1abs.eq.k2abs) then
                  if(k1abs.le.2) then
                     kf=0
                  elseif(k1abs.le.3) then
                     kf=333
                  else
                     kf=100*k1abs+10*k1abs+2*imspin+1
                  endif
               else
                  if(k1abs.gt.k2abs) then
                     kmax=k1abs
                     kmin=k2abs
                  elseif(k1abs.lt.k2abs) then
                     kmax=k2abs
                     kmin=k1abs
                  endif
                  kf=(100*kmax+10*kmin+2*imspin+1)
     1                 *isign(1,k1+k2)*(-1)**kmax
                  if(MOD(iabs(kf),10).eq.1) then
                     if(abs(sngl(xmpair)-ULMASS(iabs(kf))).gt.
     1                    abs(sngl(xmpair)-ULMASS(iabs(kf)+2))) 
     2                    kf=(iabs(kf)+2)*isign(1,kf)
                  endif
               endif
              endif
              ITYPAR(NATT)=kf
              KATT(NATT,1)=kf
            if(iabs(kf).eq.211) then
               npich=npich+1
            elseif(iabs(kf).eq.213) then
               nrhoch=nrhoch+1
            endif
           endif
 1001   CONTINUE
        if(nuudd.ne.0) then
         ppi0=float(npich/2)/float(nuudd)
         prho0=float(nrhoch/2)/float(nuudd)
      endif      
      npi0=0
      DO 1002 ISG = 1, NSG
         if(K2SGS(ISG,1).eq.-K2SGS(ISG,2)
     1        .and.iabs(K2SGS(ISG,1)).le.2.and.NJSGS(ISG).eq.2) then
            if(RANART(NSEED).le.ppi0) npi0=npi0+1
         endif
 1002 CONTINUE
      if(nuudd.gt.1) then
         call index1(MAXSTR,nuudd,xmdiag,indx)
      else
         indx(1)=1
      end if
      DO 1003 ix=1,nuudd
         iuudd=indx(ix)
         inatt=ndiag(iuudd)            
         if(ix.le.npi0) then
            kf=111
         elseif(RANART(NSEED).le.(prho0/(1-ppi0+0.00001))) then
            kf=113
         else
            if(RANART(NSEED).le.0.5) then
               kf=221
            else
               kf=223
            endif
         endif
         ITYPAR(inatt)=kf
         KATT(inatt,1)=kf
 1003 CONTINUE
      inatt=0
      if(ioscar.eq.3) then
         WRITE (85, 395) IAEVT, 3*nsmbbbar+2*nsmmeson,nsmbbbar,nsmmeson, 
     1     bimp, NELP,NINP,NELT,NINTHJ,MISS
      endif
 395  format(4I8,f10.4,5I5)
      DO 1006 ISG = 1, NSG
           if(NJSGS(ISG).ne.0) then
            inatt=inatt+1
              K1=K2SGS(ISG,1)
              k1abs=iabs(k1)
              PX1=PXSGS(ISG,1)
              PY1=PYSGS(ISG,1)
              PZ1=PZSGS(ISG,1)
              K2=K2SGS(ISG,2)
              k2abs=iabs(k2)
              PX2=PXSGS(ISG,2)
              PY2=PYSGS(ISG,2)
              PZ2=PZSGS(ISG,2)
              e1=PESGS(ISG,1)
              e2=PESGS(ISG,2)
              if(NJSGS(ISG).eq.2) then
               PXAR(inatt)=sngl(px1+px2)
               PYAR(inatt)=sngl(py1+py2)
               PZAR(inatt)=sngl(pz1+pz2)
               PATT(inatt,1)=PXAR(inatt)
               PATT(inatt,2)=PYAR(inatt)
               PATT(inatt,3)=PZAR(inatt)
               etot=e1+e2
               p1=px1+px2
               p2=py1+py2
               p3=pz1+pz2
              elseif((isoft.eq.4.or.isoft.eq.5).and.NJSGS(ISG).eq.3) 
     1              then
               PX3=PXSGS(ISG,3)
               PY3=PYSGS(ISG,3)
               PZ3=PZSGS(ISG,3)
               e3=PESGS(ISG,3)
               PXAR(inatt)=sngl(px1+px2+px3)
               PYAR(inatt)=sngl(py1+py2+py3)
               PZAR(inatt)=sngl(pz1+pz2+pz3)
               PATT(inatt,1)=PXAR(inatt)
               PATT(inatt,2)=PYAR(inatt)
               PATT(inatt,3)=PZAR(inatt)
               etot=e1+e2+e3
               p1=px1+px2+px3
               p2=py1+py2+py3
               p3=pz1+pz2+pz3
              endif
              XMAR(inatt)=ULMASS(ITYPAR(inatt))
              kf=KATT(inatt,1)
              if(kf.eq.113.or.abs(kf).eq.213.or.kf.eq.221.or.kf.eq.223
     1             .or.abs(kf).eq.313.or.abs(kf).eq.323.or.kf.eq.333
     2             .or.abs(kf).eq.1114.or.abs(kf).eq.2114
     3             .or.abs(kf).eq.2214.or.abs(kf).eq.2224) then
                 XMAR(inatt)=resmass(kf)
              endif
              PEAR(inatt)=sqrt(PXAR(inatt)**2+PYAR(inatt)**2
     1           +PZAR(inatt)**2+XMAR(inatt)**2)
              PATT(inatt,4)=PEAR(inatt)
              EATT=EATT+PEAR(inatt)
            ipartn=NJSGS(ISG)
            DO 1004 i=1,ipartn
               ftp(i)=ftsgs(isg,i)
               gxp(i)=gxsgs(isg,i)
               gyp(i)=gysgs(isg,i)
               gzp(i)=gzsgs(isg,i)
               pxp(i)=pxsgs(isg,i)
               pyp(i)=pysgs(isg,i)
               pzp(i)=pzsgs(isg,i)
               pmp(i)=pmsgs(isg,i)
               pep(i)=pesgs(isg,i)
 1004       CONTINUE
            call locldr(ipartn,drlocl)
            tau0=ARPAR1(1)
            ftavg0=ft0fom+dble(tau0)
            gxavg0=0d0
            gyavg0=0d0
            gzavg0=0d0
            DO 1005 i=1,ipartn
               gxavg0=gxavg0+gxp0(i)/ipartn
               gyavg0=gyavg0+gyp0(i)/ipartn
               gzavg0=gzavg0+gzp0(i)/ipartn
 1005       CONTINUE
            bex=p1/etot
            bey=p2/etot
            bez=p3/etot
            beta2 = bex ** 2 + bey ** 2 + bez ** 2
            gam = 1.d0 / dsqrt(1.d0 - beta2)
            if(beta2.ge.0.9999999999999d0) then
               write(6,*) '2',bex,bey,bez,beta2,gam
            endif
            call lorenz(ftavg0,gxavg0,gyavg0,gzavg0,-bex,-bey,-bez)
              GXAR(inatt)=sngl(pxnew)
              GYAR(inatt)=sngl(pynew)
              GZAR(inatt)=sngl(pznew)
              FTAR(inatt)=sngl(enenew)
              if(ioscar.eq.3) then
                 WRITE (85, 313) K2SGS(ISG,1),px1,py1,pz1,PMSGS(ISG,1),
     1                inatt,katt(inatt,1),xmar(inatt)
                 WRITE (85, 312) K2SGS(ISG,2),px2,py2,pz2,PMSGS(ISG,2),
     1                inatt,katt(inatt,1)
                 if(NJSGS(ISG).eq.3) WRITE (85, 312) K2SGS(ISG,3),
     1                px3,py3,pz3,PMSGS(ISG,3),inatt,katt(inatt,1)
              endif
 312       FORMAT(I6,4(1X,F10.3),1X,I6,1X,I6)
 313          FORMAT(I6,4(1X,F10.3),1X,I6,1X,I6,1X,F10.3)
           endif
 1006   CONTINUE
      nattzp=natt
      MSTJ(24)=mstj24
      RETURN
      END
