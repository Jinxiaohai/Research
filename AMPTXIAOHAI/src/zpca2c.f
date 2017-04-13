        subroutine zpca2c
        implicit double precision (a-h, o-z)
        character*8 code, versn
        character*4 reffra
        integer aproj, zproj, atarg, ztarg, event
        parameter (MAXPTN=400001)
        common /para1/ mul
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
        SAVE   
        data nff/0/
        if (nff .eq. 0) then
           write (26, 101) 'OSCAR1997A'
           write (26, 101) 'final_id_p_x'
           code = 'ZPC'
           versn = '1.0.1'
           aproj = -1
           zproj = -1
           atarg = -1
           ztarg = -1
           reffra = 'cm'
           ebeam = 0d0
           ntestp = 1
           write (26, 102) code, versn, aproj, zproj, atarg, ztarg,
     &        reffra, ebeam, ntestp
           nff = 1
           event = 1
           bimp = 0d0
           phi = 0d0
        end if
        write (26, 103) event, mul, bimp, phi
        do 99 i = 1, mul
           write (26, 104) i, ityp(i),
     &        px(i), py(i), pz(i), e(i), xmass(i),
     &        gx(i), gy(i), gz(i), ft(i)
 99         continue
         event = event + 1
 101        format (a12)
 102        format (2(a8, 2x), '(',i3, ',',i6, ')+(',i3, ',', i6, ')',
     &     2x, a4, 2x, e10.4, 2x, i8)
 103        format (i10, 2x, i10, 2x, f8.3, 2x, f8.3)
 104        format (i10, 2x, i10, 2x, 9(e12.6, 2x))
        return
        end
