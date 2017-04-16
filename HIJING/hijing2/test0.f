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
!
!
!....initialize HIJING for Au+Au collisions at c.m. energy of 200 GeV:
!      WRITE(*,*)  'random number seed'
      READ(*,*) dum,NSEED
!      WRITE(*,*) 'frame(LAB,CMS), enegy-frame'
      READ(*,*) dum,FRAME,EFRM
!      WRITE(*,*) 'Proj, Targ(A,P,PBAR)'
      READ(*,*) dum,PROJ,TARG
!      WRITE(*,*)  'A,Z of proj; A,Z of targ'
      READ(*,*) dum,IAP,IZP,IAT,IZT
!      WRITE(*,*) 'number of events'
      READ(*,*) dum,N_EVENT
!      WRITE(*,*) 'Print warning messages (0,1)'
      READ(*,*) dum,IHPR2(10)
!      WRITE(*,*) 'Keep information of decayed particles (0,1)'
      READ(*,*) dum,IHPR2(21)
!      WRITE(*,*) 'Impact parameter'
      READ(*,*) dum,BMIN,BMAX
!
!
!*** initialize HIJING
      CALL HIJSET(EFRM,FRAME,PROJ,TARG,IAP,IZP,IAT,IZT)
!        WRITE(*,*) 'Sjet=', HINT1(11),'mb','Stot=',HINT1(13),'mb'
!
!....set BMIN=0 and BMAX=0.0 for central interaction
      BMIN=0.0
      BMAX=HIPR1(35)+HIPR1(34)
	 

!....generating N_EVENT events of central AA interaction:
      DO 200 IE=1,N_EVENT
         CALL HIJING(FRAME,BMIN,BMAX)
         PRINT *, "BEGINNINGOFEVENT"
         !PRINT *, IE,NATT,EATT
         !PRINT *, EFRM,HINT1(19)
         !PRINT *, "PARTICLES"
         WRITE(*,'(i1.1,i7.6,ES14.7)') IE,NATT,EATT
         WRITE(*,'(ES9.3,ES14.7)') EFRM,HINT1(19)
         PRINT *, "PARTICLES"
         DO 300 IP=1,NATT
          
          IF(KATT(IP,1) >= 0) THEN
           WRITE(*,'(i5.5,i2.1,i2.1)') KATT(IP,1),KATT(IP,3),KATT(IP,4)
          ELSE 
           WRITE(*,'(i6.5,i2.1,i2.1)') KATT(IP,1),KATT(IP,3),KATT(IP,4)
          END IF
          
          IF (PATT(IP,1) >= 0 .AND. PATT(IP,2) >= 0) THEN
           WRITE(*,'(i7.7,ES14.7,ES14.7)') IP,PATT(IP,1),PATT(IP,2)
          ELSE IF (PATT(IP,1) < 0 .AND. PATT(IP,2) >= 0) THEN
           WRITE(*,'(i7.7,ES15.7,ES14.7)') IP,PATT(IP,1),PATT(IP,2)
          ELSE IF (PATT(IP,1) >= 0 .AND. PATT(IP,2) < 0) THEN
           WRITE(*,'(i7.7,ES14.7,ES15.7)') IP,PATT(IP,1),PATT(IP,2)
          ELSE 
           WRITE(*,'(i7.7,ES15.7,ES15.7)') IP,PATT(IP,1),PATT(IP,2)
          END IF
          
          IF (PATT(IP,3)>0) THEN
           WRITE(*,'(i7.7,ES14.7,ES14.7)') IP,PATT(IP,3),PATT(IP,4)
          ELSE 
           WRITE(*,'(i7.7,ES15.7,ES14.7)') IP,PATT(IP,3),PATT(IP,4)
          END IF
          !WRITE(*,'(i7.7,ES15.8,ES15.8)') IP,PATT(IP,3),PATT(IP,4)
300      continue
200   continue
      WRITE(*,*) "END"
      STOP
      END
      