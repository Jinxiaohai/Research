        double precision function ftime1(iseed)
        implicit double precision (a-h, o-z)
        external ran1
        parameter (hbarc = 0.197327054d0)
        common /par1/ formt
        SAVE   
        aa = hbarc / formt
        ftime1 = aa * dsqrt(1d0 / ran1(iseed) - 1d0)
        return
        end
