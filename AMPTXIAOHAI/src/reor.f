        subroutine reor(t, tmin, j, last0)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para5/ iconfg, iordsc
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
        SAVE   
        icels0 = icels(j)
        i1 = icels0 / 10000
        i2 = (icels0 - i1 * 10000) / 100
        i3 = icels0 - i1 * 10000 - i2 * 100
        call wallc(j, i1, i2, i3, t, tmin1)
        if (tmin .le. tmin1) then
           nc = last0
        else
           tmin = tmin1
           nc = 0
        end if
        if (iconfg .eq. 3 .or. iconfg .eq. 5) then
           call chcell(j, i1, i2, i3, last0, t, tmin, nc)
        else
           if (i1 .eq. 11 .and. i2 .eq. 11 .and. i3 .eq. 11) then
              call chout(j, last0, t, tmin, nc)
           else
              if (iconfg .eq. 1) then
                 call chin1(j, i1, i2, i3, last0, t, tmin, nc)
              else if (iconfg .eq. 2) then
                 call chin2(j, i1, i2, i3, last0, t, tmin, nc)
              else if (iconfg .eq. 4) then
                 call chin3(j, i1, i2, i3, last0, t, tmin, nc)
              end if
           end if
        end if
        call fixtim(j, t, tmin1, tmin, nc)
        return
        end
