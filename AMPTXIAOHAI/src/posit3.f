        subroutine posit3(x, y, z)
        implicit double precision (a-h, o-z)
        external ran1
        common /ilist3/ size1, size2, size3, v1, v2, v3, size
        common /rndm3/ iseedp
        SAVE   
        iseed=iseedp
         x = 2d0 * ran1(iseed) - 1d0
        y = 2d0 * ran1(iseed) - 1d0
        z = 2d0 * ran1(iseed) - 1d0
        x = x * 5d0 * size1
        y = y * 5d0 * size2
        z = z * 5d0 * size3
        return
        end
