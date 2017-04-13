       subroutine tablem
        COMMON/TABLE/ xarray(0:1000),earray(0:1000)
      SAVE   
       ptmax=2.01
       anorm=ptdis(ptmax)
       do 10 L=0,200
       x=0.01*float(L+1)
       rr=ptdis(x)/anorm
       earray(l)=rr
       xarray(l)=x
10       continue
       RETURN
       end
