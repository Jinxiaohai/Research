      program test
      CHARACTER FRAME*8,PROJ*8,TARG*8,dum*32
      COMMON/HIPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      SAVE  /HIPARNT/
      COMMON/HIMAIN1/NATT,EATT,JATT,NT,NP,N0,N01,N10,N11
      SAVE  /HIMAIN1/
      COMMON/HIMAIN2/KATT(130000,4),PATT(130000,4)
      SAVE  /HIMAIN2/
      COMMON/RANSEED/NSEED
      SAVE  /RANSEED/
C
C       
C....initialize HIJING for Au+Au collisions at c.m. energy of 200 GeV:
      WRITE(*,*)  'random number seed'
      READ(*,*) dum,NSEED
      WRITE(*,*) 'frame(LAB,CMS), enegy-frame'
      READ(*,*)  dum,FRAME,EFRM
      WRITE(*,*) 'Proj, Targ(A,P,PBAR)'
      READ(*,*)  dum,PROJ,TARG
      WRITE(*,*)  'A,Z of proj; A,Z of targ'
      READ(*,*)  dum,IAP,IZP,IAT,IZT
      WRITE(*,*) 'number of events'
      READ(*,*) dum,N_EVENT
      WRITE(*,*) 'Print warning messages (0,1)'
      READ(*,*) dum,IHPR2(10)
      WRITE(*,*) 'bmin,bmax'
      READ(*,*) dum,BMIN,BMAX
C
C
C*** initialize HIJING
        CALL HIJSET(EFRM,FRAME,PROJ,TARG,IAP,IZP,IAT,IZT)
C        type*,'Sjet=', HINT1(11),'mb','Stot=',HINT1(13),'mb'
C
C....set BMIN=0 and BMAX=0.0 for central interaction
C      BMIN=0.0
C      BMAX=0.0
C....generating N_EVENT events of central AA interaction:
      DO 200 IE=1,N_EVENT
         CALL HIJING(FRAME,BMIN,BMAX)
         WRITE(*,*) "BEGINNINGOFEVENT"
         WRITE(*,*) IE,NATT,EATT
         DO 300 I=1,NATT
            WRITE(*,*) I," ",KATT(I,1)," ",KATT(I,3)," ",KATT(I,4)
     &    ," ",PATT(I,1)," ",PATT(I,2)," ",PATT(I,3)," ",PATT(I,4)
 300     continue
 200  continue
      STOP
      END 
