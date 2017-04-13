        subroutine inian2
        implicit double precision (a-h, o-z)
        common /para5/ iconfg, iordsc
        common /ana2/
     &     det(12), dn(12), detdy(12), detdn(12), dndy(12),
     &     det1(12), dn1(12), detdy1(12), detdn1(12), dndy1(12),
     &     det2(12), dn2(12), detdy2(12), detdn2(12), dndy2(12)
        SAVE   
        if (iconfg .le. 3) then
           do 1001 i = 1, 12
              det(i) = 0d0
              dn(i) = 0d0
              det1(i) = 0d0
              dn1(i) = 0d0
              det2(i) = 0d0
              dn2(i) = 0d0
 1001      continue
        end if
        return
        end
