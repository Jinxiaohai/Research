      SUBROUTINE ARTOUT(NEVNT)
      PARAMETER (MAXSTR=150001, MAXR=1)
      PARAMETER (YMT1 = -1.0, YMT2 = 1.0)
      PARAMETER (BMT = 0.05, BY = 0.4)
      COMMON /RUN/ NUM
      COMMON /ARPRC1/ITYP1(MAXSTR, MAXR),
     &     GX1(MAXSTR, MAXR), GY1(MAXSTR, MAXR), GZ1(MAXSTR, MAXR), 
     &     FT1(MAXSTR, MAXR),
     &     PX1(MAXSTR, MAXR), PY1(MAXSTR, MAXR), PZ1(MAXSTR, MAXR),
     &     EE1(MAXSTR, MAXR), XM1(MAXSTR, MAXR)
      COMMON /ARANA1/
     &     dy1ntb(50), dy1ntp(50), DY1HM(50), 
     &     DY1KP(50), DY1KM(50), DY1K0S(50),
     &     DY1LA(50), DY1LB(50), DY1PHI(50),
     &     dm1pip(50), dm1pim(50), DMT1PR(50),
     &     DMT1PB(50), DMT1KP(50), dm1km(50),
     &     dm1k0s(50), DMT1LA(50), DMT1LB(50),
     &     dy1msn(50), DY1PIP(50), DY1PIM(50), 
     &     DY1PI0(50), DY1PR(50), DY1PB(50)
     &     ,DY1NEG(50), DY1CH(50), DE1NEG(50), DE1CH(50)
      COMMON /ARANA2/
     &     dy2ntb(50), dy2ntp(50), DY2HM(50), 
     &     DY2KP(50), DY2KM(50), DY2K0S(50),
     &     DY2LA(50), DY2LB(50), DY2PHI(50),
     &     dm2pip(50), dm2pim(50), DMT2PR(50),
     &     DMT2PB(50), DMT2KP(50), dm2km(50),
     &     dm2k0s(50), DMT2LA(50), DMT2LB(50),
     &     dy2msn(50), DY2PIP(50), DY2PIM(50), 
     &     DY2PI0(50), DY2PR(50), DY2PB(50)
     &     ,DY2NEG(50), DY2CH(50), DE2NEG(50), DE2CH(50)
      SAVE   
      OPEN (30, FILE = 'ana/dndy_netb.dat', STATUS = 'UNKNOWN')
      OPEN (31, FILE = 'ana/dndy_netp.dat', STATUS = 'UNKNOWN')
      OPEN (32, FILE = 'ana/dndy_nb.dat', STATUS = 'UNKNOWN')
      OPEN (33, FILE = 'ana/dndy_neg.dat', STATUS = 'UNKNOWN')
      OPEN (34, FILE = 'ana/dndy_ch.dat', STATUS = 'UNKNOWN')
      OPEN (35, FILE = 'ana/dnde_neg.dat', STATUS = 'UNKNOWN')
      OPEN (36, FILE = 'ana/dnde_ch.dat', STATUS = 'UNKNOWN')
      OPEN (37, FILE = 'ana/dndy_kp.dat', STATUS = 'UNKNOWN')
      OPEN (38, FILE = 'ana/dndy_km.dat', STATUS = 'UNKNOWN')
      OPEN (39, FILE = 'ana/dndy_k0l.dat', STATUS = 'UNKNOWN')
      OPEN (40, FILE = 'ana/dndy_la.dat', STATUS = 'UNKNOWN')
      OPEN (41, FILE = 'ana/dndy_lb.dat', STATUS = 'UNKNOWN')
      OPEN (42, FILE = 'ana/dndy_phi.dat', STATUS = 'UNKNOWN')
      OPEN (43, FILE = 'ana/dndy_meson.dat', STATUS = 'UNKNOWN')
      OPEN (44, FILE = 'ana/dndy_pip.dat', STATUS = 'UNKNOWN')
      OPEN (45, FILE = 'ana/dndy_pim.dat', STATUS = 'UNKNOWN')
      OPEN (46, FILE = 'ana/dndy_pi0.dat', STATUS = 'UNKNOWN')
      OPEN (47, FILE = 'ana/dndy_pr.dat', STATUS = 'UNKNOWN')
      OPEN (48, FILE = 'ana/dndy_pb.dat', STATUS = 'UNKNOWN')
      OPEN (50, FILE = 'ana/dndmtdy_pip.dat', STATUS = 'UNKNOWN')
      OPEN (51, FILE = 'ana/dndmtdy_0_1_pim.dat', STATUS = 'UNKNOWN')
      OPEN (52, FILE = 'ana/dndmtdy_pr.dat', STATUS = 'UNKNOWN')
      OPEN (53, FILE = 'ana/dndmtdy_pb.dat', STATUS = 'UNKNOWN')
      OPEN (54, FILE = 'ana/dndmtdy_kp.dat', STATUS = 'UNKNOWN')
      OPEN (55, FILE = 'ana/dndmtdy_0_5_km.dat', STATUS = 'UNKNOWN')
      OPEN (56, FILE = 'ana/dndmtdy_k0s.dat', STATUS = 'UNKNOWN')
      OPEN (57, FILE = 'ana/dndmtdy_la.dat', STATUS = 'UNKNOWN')
      OPEN (58, FILE = 'ana/dndmtdy_lb.dat', STATUS = 'UNKNOWN')
      SCALE1 = 1. / REAL(NEVNT * NUM) / BY
      SCALE2 = 1. / REAL(NEVNT * NUM) / BMT / (YMT2 - YMT1)
      DO 1001 I = 1, 50
         ymid=-10.+BY * (I - 0.5)
         WRITE (30, 333) ymid, SCALE1 * dy1ntb(I)
         WRITE (31, 333) ymid, SCALE1 * dy1ntp(I)
         WRITE (32, 333) ymid, SCALE1 * DY1HM(I)
         WRITE (37, 333) ymid, SCALE1 * DY1KP(I)
         WRITE (38, 333) ymid, SCALE1 * DY1KM(I)
         WRITE (39, 333) ymid, SCALE1 * DY1K0S(I)
         WRITE (40, 333) ymid, SCALE1 * DY1LA(I)
         WRITE (41, 333) ymid, SCALE1 * DY1LB(I)
         WRITE (42, 333) ymid, SCALE1 * DY1PHI(I)
         WRITE (33, 333) ymid, SCALE1 * DY1NEG(I)
         WRITE (34, 333) ymid, SCALE1 * DY1CH(I)
         WRITE (35, 333) ymid, SCALE1 * DE1NEG(I)
         WRITE (36, 333) ymid, SCALE1 * DE1CH(I)
         WRITE (43, 333) ymid, SCALE1 * dy1msn(I)
         WRITE (44, 333) ymid, SCALE1 * DY1PIP(I)
         WRITE (45, 333) ymid, SCALE1 * DY1PIM(I)
         WRITE (46, 333) ymid, SCALE1 * DY1PI0(I)
         WRITE (47, 333) ymid, SCALE1 * DY1PR(I)
         WRITE (48, 333) ymid, SCALE1 * DY1PB(I)
         IF (dm1pip(I) .NE. 0.0) THEN
            WRITE (50, 333) BMT * (I - 0.5), SCALE2 * dm1pip(I)
         END IF
         IF (dm1pim(I) .NE. 0.0) THEN
            WRITE (51, 333) BMT * (I - 0.5), SCALE2 * 0.1 * 
     &         dm1pim(I)
         END IF
         IF (DMT1PR(I) .NE. 0.0) THEN
            WRITE (52, 333) BMT * (I - 0.5), SCALE2 * DMT1PR(I)
         END IF
         IF (DMT1PB(I) .NE. 0.0) THEN
            WRITE (53, 333) BMT * (I - 0.5), SCALE2 * DMT1PB(I)
         END IF
         IF (DMT1KP(I) .NE. 0.0) THEN
            WRITE (54, 333) BMT * (I - 0.5), SCALE2 * DMT1KP(I)
         END IF
         IF (dm1km(I) .NE. 0.0) THEN
            WRITE (55, 333) BMT * (I - 0.5), SCALE2 * 0.5 *
     &         dm1km(I)
         END IF
         IF (dm1k0s(I) .NE. 0.0) THEN
            WRITE (56, 333) BMT * (I - 0.5), SCALE2 * dm1k0s(I)
         END IF
         IF (DMT1LA(I) .NE. 0.0) THEN
            WRITE (57, 333) BMT * (I - 0.5), SCALE2 * DMT1LA(I)
         END IF
         IF (DMT1LB(I) .NE. 0.0) THEN
            WRITE (58, 333) BMT * (I - 0.5), SCALE2 * DMT1LB(I)
         END IF
 1001 CONTINUE
      DO 1002 I = 30, 48
         WRITE (I, *) 'after hadron evolution'
 1002    CONTINUE
      DO 1003 I = 50, 58
         WRITE (I, *) 'after hadron evolution'
 1003 CONTINUE
      DO 1004 I = 1, 50
         ymid=-10.+BY * (I - 0.5)
         WRITE (30, 333) ymid, SCALE1 * dy2ntb(I)
         WRITE (31, 333) ymid, SCALE1 * dy2ntp(I)
         WRITE (32, 333) ymid, SCALE1 * DY2HM(I)
         WRITE (37, 333) ymid, SCALE1 * DY2KP(I)
         WRITE (38, 333) ymid, SCALE1 * DY2KM(I)
         WRITE (39, 333) ymid, SCALE1 * DY2K0S(I)
         WRITE (40, 333) ymid, SCALE1 * DY2LA(I)
         WRITE (41, 333) ymid, SCALE1 * DY2LB(I)
         WRITE (42, 333) ymid, SCALE1 * DY2PHI(I)
         WRITE (33, 333) ymid, SCALE1 * DY2NEG(I)
         WRITE (34, 333) ymid, SCALE1 * DY2CH(I)
         WRITE (35, 333) ymid, SCALE1 * DE2NEG(I)
         WRITE (36, 333) ymid, SCALE1 * DE2CH(I)
         WRITE (43, 333) ymid, SCALE1 * dy2msn(I)
         WRITE (44, 333) ymid, SCALE1 * DY2PIP(I)
         WRITE (45, 333) ymid, SCALE1 * DY2PIM(I)
         WRITE (46, 333) ymid, SCALE1 * DY2PI0(I)
         WRITE (47, 333) ymid, SCALE1 * DY2PR(I)
         WRITE (48, 333) ymid, SCALE1 * DY2PB(I)
         IF (dm2pip(I) .NE. 0.0) THEN
            WRITE (50, 333) BMT * (I - 0.5), SCALE2 * dm2pip(I)
         END IF
         IF (dm2pim(I) .NE. 0.0) THEN
            WRITE (51, 333) BMT * (I - 0.5), SCALE2 * 0.1 * 
     &         dm2pim(I)
         END IF
         IF (DMT2PR(I) .NE. 0.0) THEN
            WRITE (52, 333) BMT * (I - 0.5), SCALE2 * DMT2PR(I)
         END IF
         IF (DMT2PB(I) .NE. 0.0) THEN
            WRITE (53, 333) BMT * (I - 0.5), SCALE2 * DMT2PB(I)
         END IF
         IF (DMT2KP(I) .NE. 0.0) THEN
            WRITE (54, 333) BMT * (I - 0.5), SCALE2 * DMT2KP(I)
         END IF
         IF (dm2km(I) .NE. 0.0) THEN
            WRITE (55, 333) BMT * (I - 0.5), SCALE2 * 0.5 * 
     &         dm2km(I)
         END IF
         IF (dm2k0s(I) .NE. 0.0) THEN
            WRITE (56, 333) BMT * (I - 0.5), SCALE2 * dm2k0s(I)
         END IF
         IF (DMT2LA(I) .NE. 0.0) THEN
            WRITE (57, 333) BMT * (I - 0.5), SCALE2 * DMT2LA(I)
         END IF
         IF (DMT2LB(I) .NE. 0.0) THEN
            WRITE (58, 333) BMT * (I - 0.5), SCALE2 * DMT2LB(I)
         END IF
 1004 CONTINUE
 333  format(2(f12.5,1x))
      RETURN
      END
