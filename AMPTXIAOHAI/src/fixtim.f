        subroutine fixtim(l, t, tmin1, tmin, nc)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
        SAVE   
        k = nc
        if (tmin .lt. tmin1) then
           ot(l) = tmin
           if (ct(l) .lt. tmin1) then
              icsta(l) = 0
           else
              icsta(l) = icsta(l) + 10
           end if
           next(l) = k
        else if (tmin .eq. tmin1) then
           ot(l) = tmin
           if (nc .eq. 0) then
              next(l) = 0
           else
              icsta(l) = icsta(l) + 10
              next(l) = k
           end if
        else
           ot(l) = tmin1
           next(l) = 0
        end if
        return
        end
