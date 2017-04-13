        subroutine minijet_out(BB,phiRP)
        PARAMETER (MAXSTR=150001)
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
        COMMON/hjcrdn/YP(3,300),YT(3,300)
        COMMON/HJJET1/NPJ(300),KFPJ(300,500),PJPX(300,500),
     &                PJPY(300,500),PJPZ(300,500),PJPE(300,500),
     &                PJPM(300,500),NTJ(300),KFTJ(300,500),
     &                PJTX(300,500),PJTY(300,500),PJTZ(300,500),
     &                PJTE(300,500),PJTM(300,500)
        COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &       K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &       PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
        COMMON /AREVT/ IAEVT, IARUN, MISS
        common /para7/ ioscar,nsmbbbar,nsmmeson
        common/phidcy/iphidcy,pttrig,ntrig,maxmiss,ipi0dcy
        SAVE
        ntrig=0
        do I = 1, IHNT2(1)
           do J = 1, NPJ(I)
              pt=sqrt(PJPX(I,J)**2+PJPY(I,J)**2)
              if(pt.ge.pttrig) ntrig=ntrig+1
           enddo
        enddo
        do I = 1, IHNT2(3)
           do J = 1, NTJ(I)
              pt=sqrt(PJTX(I,J)**2+PJTY(I,J)**2)
              if(pt.ge.pttrig) ntrig=ntrig+1
           enddo
        enddo
        do I = 1, NSG
           do J = 1, NJSG(I)
              pt=sqrt(PXSG(I,J)**2+PYSG(I,J)**2)
              if(pt.ge.pttrig) ntrig=ntrig+1
           enddo
        enddo
        if(ntrig.eq.0) return
        if(ioscar.eq.3) write(96,*) IAEVT,MISS,IHNT2(1),IHNT2(3)
        DO 1008 I = 1, IHNT2(1)
           DO 1007 J = 1, NPJ(I)
              ityp=KFPJ(I,J)
              gx=YP(1,I)+0.5*BB*cos(phiRP)
              gy=YP(2,I)+0.5*BB*sin(phiRP)
              gz=0.
              ft=0.
              px=PJPX(I,J)
              py=PJPY(I,J)
              pz=PJPZ(I,J)
              xmass=PJPM(I,J)
              if(ioscar.eq.3) then
                 if(amax1(abs(gx),abs(gy),
     1                abs(gz),abs(ft)).lt.9999) then
                    write(96,200) ityp,px,py,pz,xmass,gx,gy,gz,ft,1
                 else
                    write(96,201) ityp,px,py,pz,xmass,gx,gy,gz,ft,1
                 endif
              endif
 1007      CONTINUE
 1008   CONTINUE
        DO 1010 I = 1, IHNT2(3)
           DO 1009 J = 1, NTJ(I)
              ityp=KFTJ(I,J)
              gx=YT(1,I)-0.5*BB*cos(phiRP)
              gy=YT(2,I)-0.5*BB*sin(phiRP)
              gz=0.
              ft=0.
              px=PJTX(I,J)
              py=PJTY(I,J)
              pz=PJTZ(I,J)
              xmass=PJTM(I,J)
              if(ioscar.eq.3) then
                 if(amax1(abs(gx),abs(gy),
     1                abs(gz),abs(ft)).lt.9999) then
                    write(96,200) ityp,px,py,pz,xmass,gx,gy,gz,ft,2
                 else
                    write(96,201) ityp,px,py,pz,xmass,gx,gy,gz,ft,2
                 endif
              endif
 1009      CONTINUE
 1010   CONTINUE
        DO 1012 I = 1, NSG
           DO 1011 J = 1, NJSG(I)
              ityp=K2SG(I,J)
              gx=0.5*(YP(1,IASG(I,1))+YT(1,IASG(I,2)))
              gy=0.5*(YP(2,IASG(I,1))+YT(2,IASG(I,2)))
              gz=0.
              ft=0.
              px=PXSG(I,J)
              py=PYSG(I,J)
              pz=PZSG(I,J)
              xmass=PMSG(I,J)
              if(ioscar.eq.3) then
                 if(amax1(abs(gx),abs(gy),
     1                abs(gz),abs(ft)).lt.9999) then
                    write(96,200) ityp,px,py,pz,xmass,gx,gy,gz,ft,3
                 else
                    write(96,201) ityp,px,py,pz,xmass,gx,gy,gz,ft,3
                 endif
              endif
 1011      CONTINUE
 1012   CONTINUE
 200  format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,2(1x,f8.2),2(2x,f2.0),2x,I2)
 201  format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,2(1x,e8.2),2(2x,f2.0),2x,I2)
        return
        end
