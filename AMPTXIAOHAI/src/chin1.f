        subroutine chin1(l, i1, i2, i3, last0, t, tmin, nc)
        implicit double precision (a-h, o-z)
        SAVE   
        itest = 0
        do 1003 i = i1 - 1, i1 + 1
           do 1002 j = i2 - 1, i2 + 1
              do 1001 k = i3 - 1, i3 + 1
                 if (i .ge. 1 .and. i .le. 10
     &              .and. j .ge. 1 .and. j .le. 10
     &              .and. k .ge. 1 .and. k .le. 10) then
                    call chcell(l, i, j, k, last0, t, tmin, nc)
                 else if (itest .eq. 0) then
                    m1 = 11
                    m2 = 11
                    m3 = 11
                    call chcell(l, m1, m2, m3, last0, t, tmin, nc)
                    itest = 1
                 end if   
 1001         continue
 1002      continue
 1003   continue
        return
        end
