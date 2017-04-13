        subroutine hoscar
        parameter (MAXSTR=150001,AMN=0.939457,AMP=0.93828)
        character*8 code, reffra, FRAME
        character*25 amptvn
        common/snn/efrm,npart1,npart2,epsiPz,epsiPt,PZPROJ,PZTARG
        common /lastt/itimeh,bimp 
        COMMON/hbt/lblast(MAXSTR),xlast(4,MAXSTR),plast(4,MAXSTR),nlast
        common/oscar1/iap,izp,iat,izt
        common/oscar2/FRAME,amptvn
        SAVE   
        data nff/0/
        if(nff.eq.0) then
           write (19, 101) 'OSCAR1997A'
           write (19, 111) 'final_id_p_x'
           code = 'AMPT'
           if(FRAME.eq.'CMS') then
              reffra = 'nncm'
              xmp=(amp*izp+amn*(iap-izp))/iap
              xmt=(amp*izt+amn*(iat-izt))/iat
              ebeam=(efrm**2-xmp**2-xmt**2)/2./xmt
           elseif(FRAME.eq.'LAB') then
              reffra = 'lab'
              ebeam=efrm
           else
              reffra = 'unknown'
              ebeam=0.
           endif
           ntestp = 1
           write (19, 102) code, amptvn, iap, izp, iat, izt,
     &        reffra, ebeam, ntestp
           nff = 1
           ievent = 1
           phi = 0.
           if(FRAME.eq.'CMS') write(19,112) efrm
        endif
        write (19, 103) ievent, nlast, bimp, phi
        do 99 i = 1, nlast
           ene=sqrt(plast(1,i)**2+plast(2,i)**2+plast(3,i)**2
     1          +plast(4,i)**2)
           write (19, 104) i, INVFLV(lblast(i)), plast(1,i),
     1          plast(2,i),plast(3,i),ene,plast(4,i),
     2          xlast(1,i),xlast(2,i),xlast(3,i),xlast(4,i)
 99     continue
        ievent = ievent + 1
 101        format (a10)
 111        format (a12)
 102        format (a4,1x,a20,1x,'(', i3, ',', i3, ')+(', i3, ',', 
     &           i3, ')', 2x, a4, 2x, e10.4, 2x, i8)
 103        format (i10, 2x, i10, 2x, f8.3, 2x, f8.3)
 104        format (i10, 2x, i10, 2x, 9(e12.6, 2x))
 112        format ('# Center-of-mass energy/nucleon-pair is',
     & f12.3,'GeV')
        return
        end
