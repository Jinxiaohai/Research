        subroutine zpca1c(p0, p1, p2, p3, ian)
        implicit double precision (a-h, o-z)
        common /ana3/ em(4, 4, 12)
        dimension en(4)
        SAVE   
        en(1) = p0
        en(2) = p1
        en(3) = p2
        en(4) = p3
        do 1002 i = 1, 4
           do 1001 j = 1, 4
              em(i, j, ian) = em(i, j, ian) + en(i) * en(j) / p0
 1001      continue
 1002   continue
        return
        end
