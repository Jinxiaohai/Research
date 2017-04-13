        integer function integ(x)
        implicit double precision (a-h, o-z)
        SAVE   
        if (x .lt. 0d0) then
           integ = int(x - 1d0)
        else
           integ = int( x )
        end if
        return
        end
