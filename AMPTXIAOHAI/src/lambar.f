         subroutine lambar(i1,i2,srt,siglab)
        PARAMETER (MAXSTR=150001)
        COMMON /AA/ R(3,MAXSTR)
        COMMON /BB/ P(3,MAXSTR)
        COMMON /CC/ E(MAXSTR)
        COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
      SAVE   
          siglab=1.e-06
        if( iabs(lb(i1)).ge.14.and.iabs(lb(i1)).le.17 )then
          eml = e(i1)
          emb = e(i2)
         else
          eml = e(i2)
          emb = e(i1)
        endif
       pthr = srt**2-eml**2-emb**2
        if(pthr .gt. 0.)then
       plab2=(pthr/2./emb)**2-eml**2
       if(plab2.gt.0)then
         plab=sqrt(plab2)
         siglab=12. + 0.43/(plab**3.3)
       if(siglab.gt.200.)siglab=200.
       endif
       endif
         return
      END
