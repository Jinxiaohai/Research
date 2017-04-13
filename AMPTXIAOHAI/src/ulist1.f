        subroutine ulist1(l, t)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para5/ iconfg, iordsc
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
        SAVE   
        icels0 = icels(l)
        i1 = icels0 / 10000
        i2 = (icels0 - i1 * 10000) / 100
        i3 = icels0 - i1 * 10000 - i2 * 100
        k = mod(icsta(l), 10)
        call wallc(l, i1, i2, i3, t, tmin1)
        tmin = tmin1
        nc = 0
        if (i1 .eq. 11 .and. i2 .eq. 11 .and. i3 .eq. 11) then
           call chkout(l, t, tmin, nc)
        else
           if (iconfg .eq. 1) then
              call chkin1(l, i1, i2, i3, t, tmin, nc)
           else if (iconfg .eq. 2) then
              call chkin2(l, i1, i2, i3, t, tmin, nc)
           else if (iconfg .eq. 4) then
              call chkin3(l, i1, i2, i3, t, tmin, nc)
           else if (iconfg .eq. 3 .or. iconfg .eq. 5) then
              call chkcel(l, i1, i2, i3, t, tmin, nc)
           end if
        end if
        call fixtim(l, t, tmin1, tmin, nc)
        return
        end
