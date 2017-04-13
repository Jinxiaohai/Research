      real function reab(i1,i2,srt,ictrl)
      PARAMETER (MAXSTR=150001,MAXR=1,PI=3.1415926)
      parameter      (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
      PARAMETER      (AKA=0.498,ALA=1.1157,ASA=1.1974)
      parameter      (amn=0.938,ap1=0.14,arho=0.77,aomega=0.782)
       parameter       (maxx=20,maxz=24)
      COMMON   /AA/  R(3,MAXSTR)
      COMMON   /BB/  P(3,MAXSTR)
      COMMON   /CC/  E(MAXSTR)
      COMMON  /DD/      RHO(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHOP(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHON(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ)
      COMMON  /EE/      ID(MAXSTR),LB(MAXSTR)
      SAVE   
       LB1=LB(I1)
       LB2=LB(I2)
       reab=0
       if(ictrl.eq.1.and.srt.le.(amn+2.*ap1+0.02))return
       if(ictrl.eq.3.and.srt.le.(amn+ap1+aomega+0.02))return
       pin2=((srt**2+ap1**2-amn**2)/(2.*srt))**2-ap1**2
       if(pin2.le.0)return
       if(ictrl.eq.1)then
       if(e(i1).gt.1)then 
       ed=e(i1)       
       else
       ed=e(i2)
       endif       
       pout2=((srt**2+ap1**2-ed**2)/(2.*srt))**2-ap1**2
       if(pout2.le.0)return
       xpro=twopi(srt)/10.
       factor=1/3.
       if( ((lb1.eq.8.and.lb2.eq.5).or.
     &    (lb1.eq.5.and.lb2.eq.8))
     &        .OR.((lb1.eq.-8.and.lb2.eq.3).or.
     &    (lb1.eq.3.and.lb2.eq.-8)) )factor=1/4.
       if((iabs(lb1).ge.10.and.iabs(lb1).le.13).
     &  or.(iabs(lb2).ge.10.and.iabs(lb2).le.13))factor=1.
       reab=factor*pin2/pout2*xpro
       return
       endif
       if(ictrl.eq.2)then
       if(lb(i2).ge.25)then 
       ed=e(i1)
       arho1=e(i2)       
       else
       ed=e(i2)
       arho1=e(i1)
       endif       
       if(srt.le.(amn+ap1+arho1+0.02))return
       pout2=((srt**2+arho1**2-ed**2)/(2.*srt))**2-arho1**2
       if(pout2.le.0)return
       xpro=threpi(srt)/10.
       factor=1/3.
       if( ((lb1.eq.8.and.lb2.eq.27).or.
     &       (lb1.eq.27.and.lb2.eq.8))
     & .OR. ((lb1.eq.-8.and.lb2.eq.25).or.
     &       (lb1.eq.25.and.lb2.eq.-8)) )factor=1/4.
       if((iabs(lb1).ge.10.and.iabs(lb1).le.13).
     &  or.(iabs(lb2).ge.10.and.iabs(lb2).le.13))factor=1.
       reab=factor*pin2/pout2*xpro
       return
       endif
       if(ictrl.eq.3)then
       if(e(i1).gt.1)ed=e(i1)       
       if(e(i2).gt.1)ed=e(i2)       
       pout2=((srt**2+aomega**2-ed**2)/(2.*srt))**2-aomega**2
       if(pout2.le.0)return
       xpro=fourpi(srt)/10.
       factor=1/6.
       if((iabs(lb1).ge.10.and.iabs(lb1).le.13).
     &  or.(iabs(lb2).ge.10.and.iabs(lb2).le.13))factor=1./3.
       reab=factor*pin2/pout2*xpro
       endif
      return
        END
