      SUBROUTINE ZPCMN
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para3/ nsevt, nevnt, nsbrun, ievt, isbrun
        SAVE   
        do 1000 i = 1, nevnt
           ievt = i
           call inievt
           do 2000 j = 1, nsbrun
              isbrun = j
              call inirun
 3000         continue
              call zpcrun(*4000)
              call zpca1
              goto 3000
 4000         continue
              call zpca2
 2000      continue
 1000   continue
        call zpcou
        CALL zpstrg
        RETURN
        end
