        subroutine inievt
        implicit double precision (a-h, o-z)
        COMMON /para1/ mul
        common /para4/ iftflg, ireflg, igeflg, ibstfg
        SAVE   
        if (ireflg .eq. 0) call readi
        if (igeflg .ne. 0) call genei
        if (ibstfg .ne. 0) call boosti
        return
        end
