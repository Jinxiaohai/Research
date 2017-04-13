        subroutine chin2(l, i1, i2, i3, last0, t, tmin, nc)
        implicit double precision (a-h, o-z)
        SAVE   
        itest = 0
        do 1003 i = i1 - 1, i1 + 1
           do 1002 j = i2 - 1, i2 + 1
              do 1001 k = i3 - 1, i3 + 1
                 ia = i
                 ib = j
                 ic = k
                 if (k .ge. 1 .and. k .le. 10) then
                    if (i .eq. 0) ia = 10
                    if (i .eq. 11) ia = 1
                    if (j .eq. 0) ib = 10
                    if (j .eq. 11) ib = 1
                    call chcell(l, ia, ib, ic, last0, t, tmin, nc)
                 end if
 1001         continue
 1002      continue
 1003   continue
        return
        end
