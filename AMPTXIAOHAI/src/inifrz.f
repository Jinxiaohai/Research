        subroutine inifrz
        implicit double precision  (a-h, o-z)
        PARAMETER (MAXPTN=400001)
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
        common /frzprc/ 
     &       gxfrz(MAXPTN), gyfrz(MAXPTN), gzfrz(MAXPTN), ftfrz(MAXPTN),
     &       pxfrz(MAXPTN), pyfrz(MAXPTN), pzfrz(MAXPTN), efrz(MAXPTN),
     &       xmfrz(MAXPTN), 
     &       tfrz(302), ifrz(MAXPTN), idfrz(MAXPTN), itlast
        SAVE   
        step1=0.1d0
        step2=1d0
        step3=10d0
        step4=100d0
        do 1001 it=1,101
           tfrz(it)=0d0+dble(it-1)*step1
 1001 continue
        do 1002 it=102,191
           tfrz(it)=10d0+dble(it-101)*step2
 1002   continue
        do 1003 it=192,281
           tfrz(it)=100d0+dble(it-191)*step3
 1003   continue
        do 1004 it=282,301
           tfrz(it)=1000d0+dble(it-281)*step4
 1004   continue
        tfrz(302)=tlarge
        return
        end
