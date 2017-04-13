        subroutine zpca1a(i)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para2/ xmp, xmu, alpha, rscut2, cutof2
        common /para5/ iconfg, iordsc
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
        common /prec3/gxs(MAXPTN),gys(MAXPTN),gzs(MAXPTN),fts(MAXPTN),
     &     pxs(MAXPTN), pys(MAXPTN), pzs(MAXPTN), es(MAXPTN),
     &     xmasss(MAXPTN), ityps(MAXPTN)
        common /prec5/ eta(MAXPTN), rap(MAXPTN), tau(MAXPTN)
        common /prec6/ etas(MAXPTN), raps(MAXPTN), taus(MAXPTN)
        common /ana1/ ts(12)
        SAVE   
        if (iconfg .eq. 1) then
           t1 = fts(i)
           t2 = ft(i)
           ipic = 11
        else if (iconfg .eq. 2 .or.
     &     iconfg .eq. 3) then
           t1 = taus(i)
           t2 = tau(i)
           ipic = 12
        else if (iconfg .eq. 4 .or.
     &     iconfg .eq. 5) then
           t1 = fts(i)
           t2 = ft(i)
           ipic = 12
        end if
        if (iconfg .le. 3) then
           do 1002 ian = 1, ipic
              if (t1 .le. ts(ian) .and.
     &           t2 .gt. ts(ian)) then
                 rapi = raps(i)
                 et = dsqrt(pxs(i) ** 2 + pys(i) ** 2 + xmp ** 2)
                 call zpca1b(rapi, et, ian)
              end if
 1002      continue
        else
           do 1003 ian = 1, ipic
              if (t1 .le. ts(ian) .and.
     &           t2 .gt. ts(ian)) then
                 p0 = es(i)
                 p1 = pxs(i)
                 p2 = pys(i)
                 p3 = pzs(i)
                 call zpca1c(p0, p1, p2, p3, ian)
              end if
 1003      continue
        end if
        return
        end
