        subroutine zpcou2
        implicit double precision (a-h, o-z)
        common /para3/ nsevt, nevnt, nsbrun, ievt, isbrun
        common /ilist3/ size1, size2, size3, v1, v2, v3, size
        common /ana1/ ts(12)
        common /ana3/ em(4, 4, 12)
        SAVE   
        open (28, file = 'ana4/em.dat', status = 'unknown')
        vol = 1000.d0 * size1 * size2 * size3
        ntotal = nevnt * nsbrun
        do 1002 ian = 1, 12
           write (28, *) '*** for time ', ts(ian), 'fm(s)'
           do 1001 i = 1, 4
              write (28, *) em(i, 1, ian) / vol / ntotal,
     &                        em(i, 2, ian) / vol / ntotal,
     &                        em(i, 3, ian) / vol / ntotal,
     &                        em(i, 4, ian) / vol / ntotal
 1001      continue
 1002   continue
        return
        end
