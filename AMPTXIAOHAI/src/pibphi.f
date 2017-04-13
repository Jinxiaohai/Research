      SUBROUTINE pibphi(srt,lb1,lb2,em1,em2,Xphi,xphin) 
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,AP2=0.13957,AM0=1.232,PI=3.1415926)
          PARAMETER (AKA=0.498, ALA = 1.1157, PIMASS=0.140, APHI=1.02)
        parameter (arho=0.77)
      SAVE   
       Xphi = 0.0
       xphin = 0.0
       xphid = 0.0
       if( (lb1.ge.3.and.lb1.le.5) .or.
     &     (lb2.ge.3.and.lb2.le.5) )then
       if( (iabs(lb1).ge.1.and.iabs(lb1).le.2) .or.
     &     (iabs(lb2).ge.1.and.iabs(lb2).le.2) )then
        IF (srt  .GT. (aphi+amn)) THEN
             srrt = srt - (aphi+amn)
             sig = 0.0235*srrt**(-0.519) 
          xphin=sig*1.*(srt**2-(aphi+amn)**2)*
     &           (srt**2-(aphi-amn)**2)/(srt**2-(em1+em2)**2)/
     &           (srt**2-(em1-em2)**2)
        END IF
        IF (srt  .GT. (aphi+am0)) THEN
             srrt = srt - (aphi+am0)
             sig = 0.0235*srrt**(-0.519) 
          xphid=sig*4.*(srt**2-(aphi+am0)**2)*
     &           (srt**2-(aphi-am0)**2)/(srt**2-(em1+em2)**2)/
     &           (srt**2-(em1-em2)**2)
        END IF
       else
        IF (srt  .GT. (aphi+amn)) THEN
             srrt = srt - (aphi+amn)
            if(srrt .lt. 0.7)then
             sig = 0.0119*srrt**(-0.534)
            else
             sig = 0.0130*srrt**(-0.304)
            endif      
          xphin=sig*(1./4.)*(srt**2-(aphi+amn)**2)*
     &           (srt**2-(aphi-amn)**2)/(srt**2-(em1+em2)**2)/
     &           (srt**2-(em1-em2)**2)
        END IF
        IF (srt  .GT. (aphi+am0)) THEN
             srrt = srt - (aphi+am0)
             if(srrt .lt. 0.7)then
             sig = 0.0119*srrt**(-0.534)
            else
             sig = 0.0130*srrt**(-0.304)
            endif      
          xphid=sig*1.*(srt**2-(aphi+am0)**2)*
     &           (srt**2-(aphi-am0)**2)/(srt**2-(em1+em2)**2)/
     &           (srt**2-(em1-em2)**2)
        END IF
       endif
       else
       if( (iabs(lb1).ge.1.and.iabs(lb1).le.2) .or.
     &     (iabs(lb2).ge.1.and.iabs(lb2).le.2) )then
        IF (srt  .GT. (aphi+amn)) THEN
             srrt = srt - (aphi+amn)
           if(srrt .lt. 0.7)then
             sig = 0.0166*srrt**(-0.786)
            else
             sig = 0.0189*srrt**(-0.277)
            endif
          xphin=sig*(1./3.)*(srt**2-(aphi+amn)**2)*
     &           (srt**2-(aphi-amn)**2)/(srt**2-(em1+em2)**2)/
     &           (srt**2-(em1-em2)**2)
        END IF
        IF (srt  .GT. (aphi+am0)) THEN
             srrt = srt - (aphi+am0)
           if(srrt .lt. 0.7)then
             sig = 0.0166*srrt**(-0.786)
            else
             sig = 0.0189*srrt**(-0.277)
            endif
          xphid=sig*(4./3.)*(srt**2-(aphi+am0)**2)*
     &           (srt**2-(aphi-am0)**2)/(srt**2-(em1+em2)**2)/
     &           (srt**2-(em1-em2)**2)
        END IF
       else
        IF (srt  .GT. (aphi+amn)) THEN
             srrt = srt - (aphi+amn)
            if(srrt .lt. 0.7)then
             sig = 0.0119*srrt**(-0.534)
            else
             sig = 0.0130*srrt**(-0.304)
            endif      
          xphin=sig*(1./12.)*(srt**2-(aphi+amn)**2)*
     &           (srt**2-(aphi-amn)**2)/(srt**2-(em1+em2)**2)/
     &           (srt**2-(em1-em2)**2)
        END IF
        IF (srt  .GT. (aphi+am0)) THEN
             srrt = srt - (aphi+am0)
             if(srrt .lt. 0.7)then
             sig = 0.0119*srrt**(-0.534)
            else
             sig = 0.0130*srrt**(-0.304)
            endif      
          xphid=sig*(1./3.)*(srt**2-(aphi+am0)**2)*
     &           (srt**2-(aphi-am0)**2)/(srt**2-(em1+em2)**2)/
     &           (srt**2-(em1-em2)**2)
        END IF
       endif
        END IF
         xphin = xphin/10.
         xphid = xphid/10.
         Xphi = xphin + xphid
       RETURN
        END
