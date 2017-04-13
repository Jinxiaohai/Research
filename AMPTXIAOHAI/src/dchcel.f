        subroutine dchcel(l, i, j, k, t)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
        common /ilist2/ icell, icel(10, 10, 10)
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
        SAVE   
        if (i .eq. 11 .or. j .eq. 11 .or. k .eq. 11) then
           if ( .not. (i .eq. 11 .and. j .eq. 11 .and.
     &     k .eq. 11)) stop 'cerr'
           m = icell
        else
           m = icel(i, j, k)
        end if
        if (m .eq. 0) return
        if (next(m) .eq. l) then
           tm = tlarge
           last0 = 0
           call reor(t, tm, m, last0)
        end if
        n = nic(m)
        if (n .eq. 0) return
        do 10 while(n .ne. m)
           if (next(n) .eq. l) then
              tm = tlarge
              last0 = 0
              call reor(t, tm, n, last0)
           end if
           n = nic(n)
 10        continue
        return
        end
