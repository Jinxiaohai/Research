        block data zpcbdt
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        PARAMETER (MAXSTR=150001)
        common /para1/ mul
        common /para2/ xmp, xmu, alpha, rscut2, cutof2
        common /para3/ nsevt, nevnt, nsbrun, ievt, isbrun
        common /para4/ iftflg, ireflg, igeflg, ibstfg
        common /para5/ iconfg, iordsc
        common /para6/ centy
        common /para7/ ioscar,nsmbbbar,nsmmeson
        COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &       PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &       XMASS0(MAXPTN), ITYP0(MAXPTN)
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
        common /prec3/gxs(MAXPTN),gys(MAXPTN),gzs(MAXPTN),fts(MAXPTN),
     &     pxs(MAXPTN), pys(MAXPTN), pzs(MAXPTN), es(MAXPTN),
     &     xmasss(MAXPTN), ityps(MAXPTN)
        common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
        common /prec5/ eta(MAXPTN), rap(MAXPTN), tau(MAXPTN)
        common /prec6/ etas(MAXPTN), raps(MAXPTN), taus(MAXPTN)
        common /aurec1/ jxa, jya, jza
        common /aurec2/ dgxa(MAXPTN), dgya(MAXPTN), dgza(MAXPTN)
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
        common /ilist2/ icell, icel(10,10,10)
        common /ilist3/ size1, size2, size3, v1, v2, v3, size
        common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
        common /ilist6/ t, iopern, icolln
        COMMON /ilist7/ LSTRG0(MAXPTN), LPART0(MAXPTN)
        COMMON /ilist8/ LSTRG1(MAXPTN), LPART1(MAXPTN)
        common /rndm1/ number
        common /rndm2/ iff
        common /rndm3/ iseedp
        common /ana1/ ts(12)
        common /ana2/
     &     det(12), dn(12), detdy(12), detdn(12), dndy(12),
     &     det1(12), dn1(12), detdy1(12), detdn1(12), dndy1(12),
     &     det2(12), dn2(12), detdy2(12), detdn2(12), dndy2(12)
        common /ana3/ em(4, 4, 12)
        common /ana4/ fdetdy(24), fdndy(24), fdndpt(12)
        SAVE   
        data centy/0d0/
        data number/0/
        data ts/0.11d0, 0.12d0, 0.15d0, 0.2d0, 0.3d0, 0.4d0, 0.6d0,
     &     0.8d0, 1d0, 2d0, 4d0, 6d0/
        end
