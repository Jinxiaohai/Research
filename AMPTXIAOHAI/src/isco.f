        subroutine isco(i, j, allok, tm, t1, t2)
        implicit double precision (a-h, o-z)
        common /para5/ iconfg, iordsc
        SAVE   
        logical allok
        iorder = iordsc / 10
        if (iconfg .eq. 1) then
           if (iorder .eq. 1) then
              call isco1(i, j, allok, tm, t1, t2)
           else if (iorder .eq. 2) then
              call isco2(i, j, allok, tm, t1, t2)
           else if (iorder .eq. 3) then
              call isco3(i, j, allok, tm, t1, t2)
           end if
        else if (iconfg .eq. 2 .or. iconfg .eq. 4) then
           if (iorder .eq. 1) then
              call isco4(i, j, allok, tm, t1, t2)
           else if (iorder .eq. 2) then
              call isco5(i, j, allok, tm, t1, t2)
           else if (iorder .eq. 3) then
              call isco6(i, j, allok, tm, t1, t2)
           end if
        else if (iconfg .eq. 3) then
           if (iorder .eq. 1) then
              call isco7(i, j, allok, tm, t1, t2)
           else if (iorder .eq. 2) then
              call isco8(i, j, allok, tm, t1, t2)
           else if (iorder .eq. 3) then
              call isco9(i, j, allok, tm, t1, t2)
           end if
        else if (iconfg .eq. 5) then
           if (iorder .eq. 1) then
              call isco10(i, j, allok, tm, t1, t2)
           else if (iorder .eq. 2) then
              call isco11(i, j, allok, tm, t1, t2)
           else if (iorder .eq. 3) then
              call isco12(i, j, allok, tm, t1, t2)
           end if
        end if
        return
        end
