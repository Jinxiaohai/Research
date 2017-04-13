        SUBROUTINE bbkaon(ic,SRT,PX,PY,PZ,ana,PlX,
     &  PlY,PlZ,ala,pkX,PkY,PkZ,icou1)
        common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      COMMON/RNDF77/NSEED
      SAVE   
       PI=3.1415962
       icou1=0
       aka=0.498
        ala=1.116
       if(ic.eq.2.or.ic.eq.4)ala=1.197
       ana=0.939
       if(ic.gt.2)then
       dmax=srt-aka-ala-0.02
        DM1=RMASS(DMAX,ISEED)
       ana=dm1
       endif
       t1=aka+ana+ala
       t2=ana+ala-aka
       if(srt.le.t1)then
       icou1=-1
       return
       endif
       pmax=sqrt((srt**2-t1**2)*(srt**2-t2**2))/(2.*srt)
       if(pmax.eq.0.)pmax=1.e-09
       ntry=0
1       pk=pmax*RANART(NSEED)
       ntry=ntry+1
       prob=fkaon(pk,pmax)
       if((prob.lt.RANART(NSEED)).and.(ntry.le.40))go to 1
       cs=1.-2.*RANART(NSEED)
       ss=sqrt(1.-cs**2)
       fai=2.*3.14*RANART(NSEED)
       pkx=pk*ss*cos(fai)
       pky=pk*ss*sin(fai)
       pkz=pk*cs
       ek=sqrt(aka**2+pk**2)
        eln=srt-ek
       if(eln.le.0)then
       icou1=-1
       return
       endif
       bx=-pkx/eln
       by=-pky/eln
       bz=-pkz/eln
       scheck=1.-bx**2-by**2-bz**2
       if(scheck.le.0) then
          write(99,*) 'scheck44: ', scheck
          stop
       endif
       ga=1./sqrt(scheck)
        elnc=eln/ga
       pn2=((elnc**2+ana**2-ala**2)/(2.*elnc))**2-ana**2
       if(pn2.le.0.)pn2=1.e-09
       pn=sqrt(pn2)
       csn=1.-2.*RANART(NSEED)
       ssn=sqrt(1.-csn**2)
       fain=2.*3.14*RANART(NSEED)
       px=pn*ssn*cos(fain)
       py=pn*ssn*sin(fain)
       pz=pn*csn
       en=sqrt(ana**2+pn2)
       plx=-px
       ply=-py
       plz=-pz
        PBETA  = PX*BX + PY*By+ PZ*Bz
              TRANS0  = GA * ( GA * PBETA / (GA + 1.) + En )
              Px = BX * TRANS0 + PX
              Py = BY * TRANS0 + PY
              Pz = BZ * TRANS0 + PZ
       el=sqrt(ala**2+plx**2+ply**2+plz**2)
        PBETA  = PlX*BX + PlY*By+ PlZ*Bz
              TRANS0  = GA * ( GA * PBETA / (GA + 1.) + El )
              Plx = BX * TRANS0 + PlX
              Ply = BY * TRANS0 + PlY
              Plz = BZ * TRANS0 + PlZ
             return
             end
