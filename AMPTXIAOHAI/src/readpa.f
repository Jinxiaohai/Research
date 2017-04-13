        subroutine readpa
        implicit double precision (a-h, o-z)
        external ran1
        character*50 str
        common /para2/ xmp, xmu, alpha, rscut2, cutof2
        common /para3/ nsevt, nevnt, nsbrun, ievt, isbrun
        common /para4/ iftflg, ireflg, igeflg, ibstfg
        common /para5/ iconfg, iordsc
        common /para7/ ioscar,nsmbbbar,nsmmeson
        common /ilist3/ size1, size2, size3, v1, v2, v3, size
        common /rndm1/ number
        common /rndm2/ iff
        common /rndm3/ iseedp
        SAVE   
        iseed=iseedp
        open (25, file = 'ana/zpc.res', status = 'unknown')
        if(ioscar.eq.1) then
           open (26, file = 'ana/parton.oscar', status = 'unknown')
           open (19, file = 'ana/hadron.oscar', status = 'unknown')
        endif
        xmp=0d0
        cutof2 = 4.5d0 * (alpha / xmu) ** 2
        rscut2=0.01d0
        nsevt=1
        nevnt=1
        nsbrun=1
        iftflg=0
        ireflg=1
        IF (ireflg .EQ. 0) THEN
           OPEN (27, FILE = 'zpc.inp', STATUS = 'UNKNOWN')
        END IF
        igeflg=0
        ibstfg=0
        iconfg=1
        iordsc=11
        v1=0.2d0
        v2=0.2d0
        v3=0.2d0
        size1=1.5d0
        size2=1.5d0
        size3=0.7d0
        if (size1 .eq. 0d0 .or. size2 .eq. 0d0 .or. 
     &     size3 .eq. 0d0) then
           if (size1 .ne. 0d0 .or. size2 .ne. 0d0 .or. size3 .ne. 0d0
     &        .or. v1 .ne. 0d0 .or. v2 .ne. 0d0 .or. v3 .ne. 0d0) then
              print *, 'to get rid of space division:'
              print *, 'set all sizes and vs to 0'
              stop 'chker'
           end if
        end if
        size = min(size1, size2, size3)
        iff=-1
        isedng=-iseed
        a = ran1(isedng)
        irused=2
        do 1001 i = 1, irused - 1
           iseed2=2
           a = ran1(iseed2)
 1001   continue
        if (iconfg .eq. 2 .or. iconfg .eq. 3) then
           v1 = 0d0
           v2 = 0d0
        end if
        if (iconfg .eq. 4 .or. iconfg .eq. 5) then
           v1 = 0d0
           v2 = 0d0
           v3 = 0d0
        end if
        close(5)
        return
        end
