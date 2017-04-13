        subroutine zpcou
        implicit double precision (a-h, o-z)
        common /para5/ iconfg, iordsc
        SAVE   
        if (iconfg .le. 3) then
           call zpcou1
        else
           call zpcou2
        end if
        return
        end
