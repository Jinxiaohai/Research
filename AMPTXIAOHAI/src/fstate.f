        subroutine fstate(iseed,srt,dm3,dm4,px,py,pz,iflag)
        dimension px(4), py(4), pz(4), pe(4)
        COMMON/RNDF77/NSEED
        SAVE   
        iflag=-1
        pio=3.1415926
        aka=0.498
        icount=0
        ekmax=(srt-dm3-dm4)/2.
        if(ekmax.le.aka) return
        pkmax=sqrt(ekmax**2-aka**2)
        if(dm3.le.0.0.or.dm4.le.0.0) then
           write(1,*)'error: minus mass!!!'
           return
        endif
50        continue
        icount=icount+1
        if(icount.gt.10) return
        ptkmi2=-1./4.145*alog(RANART(NSEED))
        ptkm=sqrt(ptkmi2)
3        v1=RANART(NSEED)
        v2=RANART(NSEED)
        rsq=v1**2+v2**2
        if(rsq.ge.1.0.or.rsq.le.0.) goto 3
        fac=sqrt(-2.*alog(rsq)/rsq)
        guass=v1*fac
        if(guass.ge.5.) goto 3
        xstar=guass/5.
        pzkm=pkmax*xstar
        ekm=sqrt(aka**2+pzkm**2+ptkm**2)
        if(RANART(NSEED).gt.aka/ekm) goto 50
        bbb=RANART(NSEED)
        px(3)=ptkm*cos(2.*pio*bbb)
        py(3)=ptkm*sin(2.*pio*bbb)
        if(RANART(NSEED).gt.0.5) pzkm=-1.*pzkm
        pz(3)=pzkm
        pe(3)=ekm
150        ptkpl2=-1./3.68*alog(RANART(NSEED))
        ptkp=sqrt(ptkpl2)
13        v1=RANART(NSEED)
        v2=RANART(NSEED)
        rsq=v1**2+v2**2
        if(rsq.ge.1.0.or.rsq.le.0.) goto 13
        fac=sqrt(-2.*alog(rsq)/rsq)
        guass=v1*fac
        if(guass.ge.3.25) goto 13
        xstar=guass/3.25
        pzkp=pkmax*xstar
        ekp=sqrt(aka**2+pzkp**2+ptkp**2)
        if(RANART(NSEED).gt.aka/ekp) goto 150
        bbb=RANART(NSEED)
        px(4)=ptkp*cos(2.*pio*bbb)
        py(4)=ptkp*sin(2.*pio*bbb)
        if(RANART(NSEED).gt.0.5) pzkp=-1.*pzkp
        pz(4)=pzkp
        pe(4)=ekp
        resten=srt-pe(3)-pe(4)
        restpz=-pz(3)-pz(4)
        if(resten.le.abs(restpz)) goto 50
        restms=sqrt(resten**2-restpz**2)
        if(restms.lt.(dm3+dm4)) goto 50
        ptp2=-1./2.76*alog(RANART(NSEED))
        ptp=sqrt(ptp2)
        bbb=RANART(NSEED)
        px(2)=ptp*cos(2.*pio*bbb)
        py(2)=ptp*sin(2.*pio*bbb)
        px(1)=-1.*(px(4)+px(3)+px(2))
        py(1)=-1.*(py(4)+py(3)+py(2))
        rmt3=sqrt(dm3**2+px(1)**2+py(1)**2)
        rmt4=sqrt(dm4**2+px(2)**2+py(2)**2)
        if(restms.lt.(rmt3+rmt4)) goto 50
        pzcms=sqrt((restms**2-(rmt3+rmt4)**2)*
     &             (restms**2-(rmt3-rmt4)**2))/2./restms
        if(RANART(NSEED).gt.0.5) then
           pz(1)=pzcms
           pz(2)=-pzcms
        else
           pz(1)=-pzcms
           pz(2)=pzcms
        endif
        beta=restpz/resten        
        gama=1./sqrt(1.-beta**2)
        pz(1)=pz(1)*gama + beta*gama*sqrt(rmt3**2+pz(1)**2)
        pz(2)=pz(2)*gama + beta*gama*sqrt(rmt4**2+pz(2)**2)
        pe(1)=sqrt(rmt3**2+pz(1)**2)
        pe(2)=sqrt(rmt4**2+pz(2)**2)
        iflag=1
        return
        end
