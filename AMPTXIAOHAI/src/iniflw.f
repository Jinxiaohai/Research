        subroutine iniflw(NEVNT,idd)
        PARAMETER (MAXSTR=150001, MAXR=1)
        DOUBLE PRECISION  v2i,eti,xmulti,v2mi,s2mi,xmmult,
     1       v2bi,s2bi,xbmult
        COMMON /RUN/ NUM
        COMMON /ARERC1/MULTI1(MAXR)
        COMMON /ARPRC1/ITYP1(MAXSTR, MAXR),
     &     GX1(MAXSTR, MAXR), GY1(MAXSTR, MAXR), GZ1(MAXSTR, MAXR), 
     &     FT1(MAXSTR, MAXR),
     &     PX1(MAXSTR, MAXR), PY1(MAXSTR, MAXR), PZ1(MAXSTR, MAXR),
     &     EE1(MAXSTR, MAXR), XM1(MAXSTR, MAXR)
        COMMON/iflow/v2i,eti,xmulti,v2mi,s2mi,xmmult,v2bi,s2bi,xbmult
        SAVE   
        if(idd.eq.0) then
           v2i=0d0
           eti=0d0
           xmulti=0d0
           v2mi=0d0
           s2mi=0d0
           xmmult=0d0
           v2bi=0d0
           s2bi=0d0
           xbmult=0d0
        else if(idd.eq.1) then
           do 1002 J = 1, NUM
              do 1001 I = 1, MULTI1(J)
                 ITYP = ITYP1(I, J)
                 IF (ITYP .GT. -100 .AND. ITYP .LT. 100) GOTO 100
                 xmulti=xmulti+1.d0
                 PX = PX1(I, J)
                 PY = PY1(I, J)
                 XM = XM1(I, J)
                 pt2=px**2+py**2
                 xh=gx1(I,J)
                 yh=gy1(I,J)
                 xt2=xh**2+yh**2
                 if(pt2.gt.0) v2i=v2i+dble((px**2-py**2)/pt2)
                 eti=eti+dble(SQRT(PX ** 2 + PY ** 2 + XM ** 2))
                 IF (ITYP .LT. -1000 .or. ITYP .GT. 1000) then
                    xbmult=xbmult+1.d0
                    if(pt2.gt.0) v2bi=v2bi+dble((px**2-py**2)/pt2)
                    if(xt2.gt.0) s2bi=s2bi+dble((xh**2-yh**2)/xt2)
                 else
                    xmmult=xmmult+1.d0
                    if(pt2.gt.0) v2mi=v2mi+dble((px**2-py**2)/pt2)
                    if(xt2.gt.0) s2mi=s2mi+dble((xh**2-yh**2)/xt2)
                 endif
 100                 continue
 1001         continue
 1002      continue
        else if(idd.eq.2) then
           if(xmulti.ne.0) v2i=v2i/xmulti
           eti=eti/dble(NEVNT)
           xmulti=xmulti/dble(NEVNT)
           if(xmmult.ne.0) then
              v2mi=v2mi/xmmult
              s2mi=s2mi/xmmult
           endif
           xmmult=xmmult/dble(NEVNT)
           if(xbmult.ne.0) then
              v2bi=v2bi/xbmult
              s2bi=s2bi/xbmult
           endif
           xbmult=xbmult/dble(NEVNT)
        endif
        return
        end
