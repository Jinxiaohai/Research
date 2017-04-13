        subroutine lorntz(ilo,b,pi,pj)
        dimension pi(4),pj(4),b(3)
      SAVE   
        bb=b(1)*b(1)+b(2)*b(2)+b(3)*b(3)
        deno3=sqrt(1.-bb)
        if(deno3.eq.0.)deno3=1.e-10
        gam=1./deno3
        ga=gam*gam/(gam+1.)
        if(ilo.eq.1) goto 100
        pib=pi(1)*b(1)+pi(2)*b(2)+pi(3)*b(3)
        pjb=pj(1)*b(1)+pj(2)*b(2)+pj(3)*b(3)
        do 1001 i=1,3
           pi(i)=pi(i)+b(i)*(ga*pib-gam*pi(4))
           pj(i)=pj(i)+b(i)*(ga*pjb-gam*pj(4))
 1001   continue
        pi(4)=gam*(pi(4)-pib)
        pj(4)=gam*(pj(4)-pjb)
        return
100     continue
        pib=pi(1)*b(1)+pi(2)*b(2)+pi(3)*b(3)
        pjb=pj(1)*b(1)+pj(2)*b(2)+pj(3)*b(3)
        do 1002 i=1,3
           pi(i)=pi(i)+b(i)*(ga*pib+gam*pi(4))
           pj(i)=pj(i)+b(i)*(ga*pjb+gam*pj(4))
 1002   continue
        pi(4)=gam*(pi(4)+pib)
        pj(4)=gam*(pj(4)+pjb)
        return
        end
