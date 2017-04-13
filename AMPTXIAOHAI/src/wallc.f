        subroutine wallc(i, i1, i2, i3, t, tmin)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para5/ iconfg, iordsc
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
        SAVE   
        tmin = tlarge
        if (iconfg .le. 2 .or. iconfg .eq. 4) then
           if ((i1 .ge. 1 .and. i1 .le. 10)
     &          .or. (i2 .ge. 1 .and. i2 .le. 10)
     &          .or. (i3 .ge. 1 .and. i3 .le. 10)) then
              call wallc1(i, i1, i2, i3, t, tmin)
           else
              call wallcb(i, t, tmin)              
           end if
        else if (iconfg .eq. 3 .or. iconfg .eq. 5) then
           call wallc2(i, i1, i2, i3, t, tmin)
        end if
        return
        end
