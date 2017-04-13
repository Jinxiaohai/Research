        subroutine zpca1b(rapi, et, ian)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para6/ centy
        common /ilist6/ t, iopern, icolln
        common /ana2/
     &     det(12), dn(12), detdy(12), detdn(12), dndy(12),
     &     det1(12), dn1(12), detdy1(12), detdn1(12), dndy1(12),
     &     det2(12), dn2(12), detdy2(12), detdn2(12), dndy2(12)
        SAVE   
        if (rapi .gt. centy - 0.5d0 .and. 
     &     rapi .lt. centy + 0.5d0) then
           det2(ian) = det2(ian) + et
           dn2(ian) = dn2(ian) + 1d0
           if (ian .eq. 10) then
           end if
           if (ian .eq. 11) then
           end if
           if (ian .eq. 12) then
           end if
           if (rapi .gt. centy - 0.25d0 .and. 
     &        rapi .lt. centy + 0.25d0) then
              det1(ian) = det1(ian) + et
              dn1(ian) = dn1(ian) + 1d0
              if (rapi .gt. centy - 0.1d0 .and.
     &           rapi .lt. centy + 0.1d0) then
                 det(ian) = det(ian) + et
                 dn(ian) = dn(ian) + 1d0
              end if
           end if
        end if
        return
        end
