        SUBROUTINE ppomga(SRT,ISEED,PX,PY,PZ,DM1,PNX,
     &  PNY,PNZ,DM2,PPX,PPY,PPZ,icou1)
        COMMON/TABLE/ xarray(0:1000),earray(0:1000)
      COMMON/RNDF77/NSEED
      SAVE   
        ntrym=0
       icou1=0
       pi=3.1415926
        AMN=938.925/1000.
        AMP=782./1000.
       DM1=amn
       DM2=amn
       V=0.43
       W=-0.84
       PTMAX2=(SRT**2-(DM1+DM2+AMP)**2)*
     1  (SRT**2-(DM1-AMP-DM2)**2)/4./SRT**2
       scheck=PTMAX2
       if(scheck.lt.0) then
          write(99,*) 'scheck30: ', scheck
          scheck=0.
       endif
       PTMAX=SQRT(scheck)*1./3.
7       PT=PTR(PTMAX,ISEED)
       PZMAX2=(SRT**2-(DM1+DM2+AMP)**2)*
     1  (SRT**2-(DM1-AMP-DM2)**2)/4./SRT**2-PT**2
       NTRYM=NTRYM+1
       IF((PZMAX2.LT.0.).and.ntrym.le.100)then 
       go to 7
       else
       pzmax2=1.E-09
       endif
       PZMAX=SQRT(PZMAX2)
       XMAX=2.*PZMAX/SRT
       ntryx=0
       fmax00=1.056
       x00=0.26
       if(abs(xmax).gt.0.26)then
       f00=fmax00
       else
       f00=1.+v*abs(xmax)+w*xmax**2
       endif
9       X=XMAX*(1.-2.*RANART(NSEED))
       ntryx=ntryx+1
       xratio=(1.+V*ABS(X)+W*X**2)/f00       
       IF(xratio.LT.RANART(NSEED).and.ntryx.le.50)GO TO 9       
       PZ=0.5*SRT*X
       fai=2.*pi*RANART(NSEED)
       Px=pt*cos(fai)
       Py=pt*sin(fai)
       ek=sqrt(dm1**2+PT**2+Pz**2)
        eln=srt-ek
       IF(ELN.lE.0)then
       icou1=-1
       return
       endif
       bx=-Px/eln
       by=-Py/eln
       bz=-Pz/eln
       scheck=1.-bx**2-by**2-bz**2
       if(scheck.le.0) then
          write(99,*) 'scheck31: ', scheck
          stop
       endif
       ga=1./sqrt(scheck)
       elnc=eln/ga
       pn2=((elnc**2+dm2**2-amp**2)/(2.*elnc))**2-dm2**2
       if(pn2.le.0)then
       icou1=-1
       return
       endif
       pn=sqrt(pn2)
        xptr=0.33*PN
       PNT=PTR(xptr,ISEED)
       fain=2.*pi*RANART(NSEED)
       pnx=pnT*cos(fain)
       pny=pnT*sin(fain)
       SIG=1
       IF(X.GT.0)SIG=-1
       scheck=pn**2-PNT**2
       if(scheck.lt.0) then
          write(99,*) 'scheck32: ', scheck
          scheck=0.
       endif
       pnz=SIG*SQRT(scheck)
       en=sqrt(dm2**2+pnx**2+pny**2+pnz**2)
       ppx=-pnx
       ppy=-pny
       ppz=-pnz
       ep=sqrt(amp**2+ppx**2+ppy**2+ppz**2)
        PBETA  = PnX*BX + PnY*By+ PnZ*Bz
              TRANS0  = GA * ( GA * PBETA / (GA + 1.) + En )
              Pnx = BX * TRANS0 + PnX
              Pny = BY * TRANS0 + PnY
              Pnz = BZ * TRANS0 + PnZ
             if(ep.eq.0.)ep=1.E-09
              PBETA  = PPX*BX + PPY*By+ PPZ*Bz
              TRANS0  = GA * ( GA * PBETA / (GA + 1.) + EP )
              PPx = BX * TRANS0 + PPX
              PPy = BY * TRANS0 + PPY
              PPz = BZ * TRANS0 + PPZ
       return
       end
