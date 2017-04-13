        subroutine frztm(NEVNT,idd)
        implicit double precision  (a-h, o-z)
        PARAMETER (MAXPTN=400001)
        dimension tsf(31)
        COMMON /PARA1/ MUL
        COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &       PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &       XMASS0(MAXPTN), ITYP0(MAXPTN)
        COMMON /prec2/GX5(MAXPTN),GY5(MAXPTN),GZ5(MAXPTN),FT5(MAXPTN),
     &       PX5(MAXPTN), PY5(MAXPTN), PZ5(MAXPTN), E5(MAXPTN),
     &       XMASS5(MAXPTN), ITYP5(MAXPTN)
        COMMON /frzout/ xnprod(30),etprod(30),xnfrz(30),etfrz(30),
     & dnprod(30),detpro(30),dnfrz(30),detfrz(30)
        SAVE   
        data tsf/0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
     &       1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
     &       2  , 3,   4,   5,   6,   7,   8,   9,   10,  20,  30/
        if(idd.eq.0) then
           do 1001 ii=1,30
              xnprod(ii)=0d0
              etprod(ii)=0d0
              xnfrz(ii)=0d0
              etfrz(ii)=0d0
              dnprod(ii)=0d0
              detpro(ii)=0d0
              dnfrz(ii)=0d0
              detfrz(ii)=0d0
 1001      continue
           OPEN (86, FILE = 'ana1/production.dat', STATUS = 'UNKNOWN')
           OPEN (87, FILE = 'ana1/freezeout.dat', STATUS = 'UNKNOWN')
        else if(idd.eq.1) then
           DO 100 ip = 1, MUL
              do 1002 ii=1,30
                 eth0=dSQRT(PX0(ip)**2+PY0(ip)**2+XMASS0(ip)**2)
                 eth2=dSQRT(PX5(ip)**2+PY5(ip)**2+XMASS5(ip)**2)
                 if (ft0(ip).lt.tsf(ii+1)) then
                    xnprod(ii)=xnprod(ii)+1d0
                    etprod(ii)=etprod(ii)+eth0
                    if (ft0(ip).ge.tsf(ii)) then
                       dnprod(ii)=dnprod(ii)+1d0
                       detpro(ii)=detpro(ii)+eth0
                    endif
                 endif
                 if (FT5(ip).lt.tsf(ii+1)) then
                    xnfrz(ii)=xnfrz(ii)+1d0
                    etfrz(ii)=etfrz(ii)+eth2
                    if (FT5(ip).ge.tsf(ii)) then
                       dnfrz(ii)=dnfrz(ii)+1d0
                       detfrz(ii)=detfrz(ii)+eth2
                    endif
                 endif
 1002         continue
 100           continue
        else if(idd.eq.2) then
           write (86,*) '       t,       np,       dnp/dt,      etp '//
     1 ' detp/dt'
           write (87,*) '       t,       nf,       dnf/dt,      etf '//
     1 ' detf/dt'
           do 1003 ii=1,30
              xnp=xnprod(ii)/dble(NEVNT)
              xnf=xnfrz(ii)/dble(NEVNT)
              etp=etprod(ii)/dble(NEVNT)
              etf=etfrz(ii)/dble(NEVNT)
              dxnp=dnprod(ii)/dble(NEVNT)/(tsf(ii+1)-tsf(ii))
              dxnf=dnfrz(ii)/dble(NEVNT)/(tsf(ii+1)-tsf(ii))
              detp=detpro(ii)/dble(NEVNT)/(tsf(ii+1)-tsf(ii))
              detf=detfrz(ii)/dble(NEVNT)/(tsf(ii+1)-tsf(ii))
              write (86, 200) 
     1        tsf(ii+1),xnp,dxnp,etp,detp
              write (87, 200) 
     1        tsf(ii+1),xnf,dxnf,etf,detf
 1003      continue
        endif
 200    format(2x,f9.2,4(2x,f10.2))
        return
        end
