        subroutine posit1(x, y, r0)
        implicit double precision (a-h, o-z)
        external ran1
        common /rndm3/ iseedp
        SAVE   
        iseed=iseedp
 10        x = 2d0 * ran1(iseed) - 1d0
        y = 2d0 * ran1(iseed) - 1d0
        if (x ** 2 + y ** 2 .gt. 1d0) goto 10
        x = x * r0
        y = y * r0
        return
        end
