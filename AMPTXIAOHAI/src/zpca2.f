        subroutine zpca2
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para3/ nsevt, nevnt, nsbrun, ievt, isbrun
        common /para5/ iconfg, iordsc
        common /para7/ ioscar,nsmbbbar,nsmmeson
        common /ilist6/ t, iopern, icolln
        common /rndm1/ number
        common /rndm2/ iff
        common /rndm3/ iseedp
        COMMON /AREVT/ IAEVT, IARUN, MISS
        SAVE   
        if (iconfg .le. 3) then
           call zpca2a
        else
           call zpca2b
        end if
        if (ioscar .eq. 1) then
           call zpca2c
        end if
        WRITE (25, *) ' Event ', IAEVT, ', run ', IARUN
        WRITE (25, *) '    number of operations = ', iopern
        WRITE (25, *) '    number of collisions between particles = ', 
     &       icolln
        WRITE (25, *) '    freezeout time=', t
        WRITE (25, *) '    ending at the ', number, 'th random number'
        WRITE (25, *) '    ending collision iff=', iff
        return
        end
