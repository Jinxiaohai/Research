        subroutine savrec(i)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
        common /prec3/gxs(MAXPTN),gys(MAXPTN),gzs(MAXPTN),fts(MAXPTN),
     &     pxs(MAXPTN), pys(MAXPTN), pzs(MAXPTN), es(MAXPTN),
     &     xmasss(MAXPTN), ityps(MAXPTN)
        common /prec5/ eta(MAXPTN), rap(MAXPTN), tau(MAXPTN)
        common /prec6/ etas(MAXPTN), raps(MAXPTN), taus(MAXPTN)
        SAVE   
        ityps(i) = ityp(i)
        gxs(i) = gx(i)
        gys(i) = gy(i)
        gzs(i) = gz(i)
        fts(i) = ft(i)
        pxs(i) = px(i)
        pys(i) = py(i)
        pzs(i) = pz(i)
        es(i) = e(i)
        xmasss(i) = xmass(i)
        etas(i) = eta(i)
        raps(i) = rap(i)
        taus(i) = tau(i)
        return
        end
