        subroutine energy(e, temp)
        implicit double precision (a-h, o-z)
        external ran1
        common /para2/ xmp, xmu, alpha, rscut2, cutof2
        common /rndm3/ iseedp
        SAVE   
        iseed=iseedp
 1000        continue
        e = ran1(iseed)
        e = e * ran1(iseed)
        e = e * ran1(iseed)
        if (e .le. 0d0) goto 1000
        e = - temp * log(e)
        if (ran1(iseed) .gt. 
     &     exp((e - dsqrt(e ** 2 + xmp ** 2))/temp)) then
           goto 1000
        end if
        return
        end
