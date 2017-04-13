      SUBROUTINE PYRAND 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200) 
      SAVE /PYSUBS/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2) 
      SAVE /PYINT2/ 
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)  
      SAVE /PYINT3/ 
      COMMON/PYINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3) 
      SAVE /PYINT4/ 
      COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3) 
      SAVE /PYINT5/ 
      MINT(17)=0    
      MINT(18)=0    
      VINT(143)=1.  
      VINT(144)=1.  
      IF(MSUB(95).EQ.1.OR.MINT(82).GE.2) CALL PYMULT(2) 
      ISUB=0    
  100 MINT(51)=0    
      IF(MINT(82).EQ.1.AND.(ISUB.LE.90.OR.ISUB.GT.96)) THEN 
        RSUB=XSEC(0,1)*RLU(0)   
        DO 110 I=1,200  
        IF(MSUB(I).NE.1) GOTO 110   
        ISUB=I  
        RSUB=RSUB-XSEC(I,1) 
        IF(RSUB.LE.0.) GOTO 120 
  110   CONTINUE    
  120   IF(ISUB.EQ.95) ISUB=96  
      ELSEIF(MINT(82).GE.2.AND.ISUB.EQ.0) THEN  
        RSUB=VINT(131)*RLU(0)   
        ISUB=96 
        IF(RSUB.GT.VINT(106)) ISUB=93   
        IF(RSUB.GT.VINT(106)+VINT(104)) ISUB=92 
        IF(RSUB.GT.VINT(106)+VINT(104)+VINT(103)) ISUB=91   
      ENDIF 
      IF(MINT(82).EQ.1) NGEN(0,1)=NGEN(0,1)+1   
      IF(MINT(82).EQ.1) NGEN(ISUB,1)=NGEN(ISUB,1)+1 
      MINT(1)=ISUB  
      MINT(72)=0    
      KFR1=0    
      IF(ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.3) THEN   
        KFR1=KFPR(ISUB,1)   
      ELSEIF(ISUB.GE.71.AND.ISUB.LE.77) THEN    
        KFR1=25 
      ENDIF 
      IF(KFR1.NE.0) THEN    
        TAUR1=PMAS(KFR1,1)**2/VINT(2)   
        GAMR1=PMAS(KFR1,1)*PMAS(KFR1,2)/VINT(2) 
        MINT(72)=1  
        MINT(73)=KFR1   
        VINT(73)=TAUR1  
        VINT(74)=GAMR1  
      ENDIF 
      IF(ISUB.EQ.141) THEN  
        KFR2=23 
        TAUR2=PMAS(KFR2,1)**2/VINT(2)   
        GAMR2=PMAS(KFR2,1)*PMAS(KFR2,2)/VINT(2) 
        MINT(72)=2  
        MINT(74)=KFR2   
        VINT(75)=TAUR2  
        VINT(76)=GAMR2  
      ENDIF 
      VINT(63)=0.   
      VINT(64)=0.   
      MINT(71)=0    
      VINT(71)=CKIN(3)  
      IF(MINT(82).GE.2) VINT(71)=0. 
      IF(ISET(ISUB).EQ.2.OR.ISET(ISUB).EQ.4) THEN   
        DO 130 I=1,2    
        IF(KFPR(ISUB,I).EQ.0) THEN  
        ELSEIF(MSTP(42).LE.0) THEN  
          VINT(62+I)=PMAS(KFPR(ISUB,I),1)**2    
        ELSE    
          VINT(62+I)=ULMASS(KFPR(ISUB,I))**2    
        ENDIF   
  130   CONTINUE    
        IF(MIN(VINT(63),VINT(64)).LT.CKIN(6)**2) MINT(71)=1 
        IF(MINT(71).EQ.1) VINT(71)=MAX(CKIN(3),CKIN(5)) 
      ENDIF 
      IF(ISET(ISUB).EQ.0) THEN  
        IS=INT(1.5+RLU(0))  
        VINT(63)=VINT(3)**2 
        VINT(64)=VINT(4)**2 
        IF(ISUB.EQ.92.OR.ISUB.EQ.93) VINT(62+IS)=PARP(111)**2   
        IF(ISUB.EQ.93) VINT(65-IS)=PARP(111)**2 
        SH=VINT(2)  
        SQM1=VINT(3)**2 
        SQM2=VINT(4)**2 
        SQM3=VINT(63)   
        SQM4=VINT(64)   
        SQLA12=(SH-SQM1-SQM2)**2-4.*SQM1*SQM2   
        SQLA34=(SH-SQM3-SQM4)**2-4.*SQM3*SQM4   
        THTER1=SQM1+SQM2+SQM3+SQM4-(SQM1-SQM2)*(SQM3-SQM4)/SH-SH    
        THTER2=SQRT(MAX(0.,SQLA12))*SQRT(MAX(0.,SQLA34))/SH 
        THL=0.5*(THTER1-THTER2) 
        THU=0.5*(THTER1+THTER2) 
        THM=MIN(MAX(THL,PARP(101)),THU) 
        JTMAX=0 
        IF(ISUB.EQ.92.OR.ISUB.EQ.93) JTMAX=ISUB-91  
        DO 140 JT=1,JTMAX   
        MINT(13+3*JT-IS*(2*JT-3))=1 
        SQMMIN=VINT(59+3*JT-IS*(2*JT-3))    
        SQMI=VINT(8-3*JT+IS*(2*JT-3))**2    
        SQMJ=VINT(3*JT-1-IS*(2*JT-3))**2    
        SQMF=VINT(68-3*JT+IS*(2*JT-3))  
        SQUA=0.5*SH/SQMI*((1.+(SQMI-SQMJ)/SH)*THM+SQMI-SQMF-    
     &  SQMJ**2/SH+(SQMI+SQMJ)*SQMF/SH+(SQMI-SQMJ)**2/SH**2*SQMF)   
        QUAR=SH/SQMI*(THM*(THM+SH-SQMI-SQMJ-SQMF*(1.-(SQMI-SQMJ)/SH))+  
     &  SQMI*SQMJ-SQMJ*SQMF*(1.+(SQMI-SQMJ-SQMF)/SH))   
        SQMMAX=SQUA+SQRT(MAX(0.,SQUA**2-QUAR))  
        IF(ABS(QUAR/SQUA**2).LT.1.E-06) SQMMAX=0.5*QUAR/SQUA    
        SQMMAX=MIN(SQMMAX,(VINT(1)-SQRT(SQMF))**2)  
        VINT(59+3*JT-IS*(2*JT-3))=SQMMIN*(SQMMAX/SQMMIN)**RLU(0)    
  140   CONTINUE    
        SQM3=VINT(63)   
        SQM4=VINT(64)   
        SQLA34=(SH-SQM3-SQM4)**2-4.*SQM3*SQM4   
        THTER1=SQM1+SQM2+SQM3+SQM4-(SQM1-SQM2)*(SQM3-SQM4)/SH-SH    
        THTER2=SQRT(MAX(0.,SQLA12))*SQRT(MAX(0.,SQLA34))/SH 
        THL=0.5*(THTER1-THTER2) 
        THU=0.5*(THTER1+THTER2) 
        B=VINT(121) 
        C=VINT(122) 
        IF(ISUB.EQ.92.OR.ISUB.EQ.93) THEN   
          B=0.5*B   
          C=0.5*C   
        ENDIF   
        THM=MIN(MAX(THL,PARP(101)),THU) 
        EXPTH=0.    
        THARG=B*(THM-THU)   
        IF(THARG.GT.-20.) EXPTH=EXP(THARG)  
  150   TH=THU+LOG(EXPTH+(1.-EXPTH)*RLU(0))/B   
        TH=MAX(THM,MIN(THU,TH)) 
        RATLOG=MIN((B+C*(TH+THM))*(TH-THM),(B+C*(TH+THU))*(TH-THU)) 
        IF(RATLOG.LT.LOG(RLU(0))) GOTO 150  
        VINT(21)=1. 
        VINT(22)=0. 
        VINT(23)=MIN(1.,MAX(-1.,(2.*TH-THTER1)/THTER2)) 
      ELSEIF(ISET(ISUB).GE.1.AND.ISET(ISUB).LE.4) THEN  
        CALL PYKLIM(1)  
        IF(MINT(51).NE.0) GOTO 100  
        RTAU=RLU(0) 
        MTAU=1  
        IF(RTAU.GT.COEF(ISUB,1)) MTAU=2 
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)) MTAU=3    
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)) MTAU=4   
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)+COEF(ISUB,4)) 
     &  MTAU=5  
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)+COEF(ISUB,4)+ 
     &  COEF(ISUB,5)) MTAU=6    
        CALL PYKMAP(1,MTAU,RLU(0))  
        IF(ISET(ISUB).EQ.3.OR.ISET(ISUB).EQ.4) THEN 
          CALL PYKLIM(4)    
          IF(MINT(51).NE.0) GOTO 100    
          RTAUP=RLU(0)  
          MTAUP=1   
          IF(RTAUP.GT.COEF(ISUB,15)) MTAUP=2    
          CALL PYKMAP(4,MTAUP,RLU(0))   
        ENDIF   
        CALL PYKLIM(2)  
        IF(MINT(51).NE.0) GOTO 100  
        RYST=RLU(0) 
        MYST=1  
        IF(RYST.GT.COEF(ISUB,7)) MYST=2 
        IF(RYST.GT.COEF(ISUB,7)+COEF(ISUB,8)) MYST=3    
        CALL PYKMAP(2,MYST,RLU(0))  
        CALL PYKLIM(3)  
        IF(MINT(51).NE.0) GOTO 100  
        IF(ISET(ISUB).EQ.2.OR.ISET(ISUB).EQ.4) THEN 
          RCTH=RLU(0)   
          MCTH=1    
          IF(RCTH.GT.COEF(ISUB,10)) MCTH=2  
          IF(RCTH.GT.COEF(ISUB,10)+COEF(ISUB,11)) MCTH=3    
          IF(RCTH.GT.COEF(ISUB,10)+COEF(ISUB,11)+COEF(ISUB,12)) MCTH=4  
          IF(RCTH.GT.COEF(ISUB,10)+COEF(ISUB,11)+COEF(ISUB,12)+ 
     &    COEF(ISUB,13)) MCTH=5 
          CALL PYKMAP(3,MCTH,RLU(0))    
        ENDIF   
      ELSEIF(ISET(ISUB).EQ.5) THEN  
        CALL PYMULT(3)  
        ISUB=MINT(1)    
      ENDIF 
      VINT(24)=PARU(2)*RLU(0)   
      MINT(51)=0    
      IF(ISUB.LE.90.OR.ISUB.GT.100) CALL PYKLIM(0)  
      IF(MINT(51).NE.0) GOTO 100    
      IF(MINT(82).EQ.1.AND.MSTP(141).GE.1) THEN 
        MCUT=0  
        IF(MSUB(91)+MSUB(92)+MSUB(93)+MSUB(94)+MSUB(95).EQ.0)   
     &  CALL PYKCUT(MCUT)   
        IF(MCUT.NE.0) GOTO 100  
      ENDIF 
      CALL PYSIGH(NCHN,SIGS)    
      IF(MINT(82).EQ.1.AND.ISUB.LE.90.OR.ISUB.GE.96) THEN   
        XSEC(ISUB,2)=XSEC(ISUB,2)+SIGS  
      ELSEIF(MINT(82).EQ.1) THEN    
        XSEC(ISUB,2)=XSEC(ISUB,2)+XSEC(ISUB,1)  
      ENDIF 
      IF(MINT(43).EQ.4.AND.MSTP(82).GE.3) THEN  
        VINT(153)=SIGS  
        CALL PYMULT(4)  
      ENDIF 
      VIOL=SIGS/XSEC(ISUB,1)    
      IF(VIOL.LT.RLU(0)) GOTO 100   
      IF(MSTP(123).LE.0) THEN   
        IF(VIOL.GT.1.) THEN 
          WRITE(MSTU(11),1000) VIOL,NGEN(0,3)+1 
          WRITE(MSTU(11),1100) ISUB,VINT(21),VINT(22),VINT(23),VINT(26) 
          STOP  
        ENDIF   
      ELSEIF(MSTP(123).EQ.1) THEN   
        IF(VIOL.GT.VINT(108)) THEN  
          VINT(108)=VIOL    
        ENDIF   
      ELSEIF(VIOL.GT.VINT(108)) THEN    
        VINT(108)=VIOL  
        IF(VIOL.GT.1.) THEN 
          XDIF=XSEC(ISUB,1)*(VIOL-1.)   
          XSEC(ISUB,1)=XSEC(ISUB,1)+XDIF    
          IF(MSUB(ISUB).EQ.1.AND.(ISUB.LE.90.OR.ISUB.GT.96))    
     &    XSEC(0,1)=XSEC(0,1)+XDIF  
          VINT(108)=1.  
        ENDIF   
      ENDIF 
      VINT(148)=1.  
      IF(MINT(43).EQ.4.AND.(ISUB.LE.90.OR.ISUB.GE.96).AND.MSTP(82).GE.3)    
     &THEN  
        CALL PYMULT(5)  
        IF(VINT(150).LT.RLU(0)) GOTO 100    
      ENDIF 
      IF(MINT(82).EQ.1.AND.MSUB(95).EQ.1) THEN  
        IF(ISUB.LE.90.OR.ISUB.GE.95) NGEN(95,1)=NGEN(95,1)+1    
        IF(ISUB.LE.90.OR.ISUB.GE.96) NGEN(96,2)=NGEN(96,2)+1    
      ENDIF 
      IF(ISUB.LE.90.OR.ISUB.GE.96) MINT(31)=MINT(31)+1  
      RSIGS=SIGS*RLU(0) 
      QT2=VINT(48)  
      RQQBAR=PARP(87)*(1.-(QT2/(QT2+(PARP(88)*PARP(82))**2))**2)    
      IF(ISUB.NE.95.AND.(ISUB.NE.96.OR.MSTP(82).LE.1.OR.    
     &RLU(0).GT.RQQBAR)) THEN   
        DO 190 ICHN=1,NCHN  
        KFL1=ISIG(ICHN,1)   
        KFL2=ISIG(ICHN,2)   
        MINT(2)=ISIG(ICHN,3)    
        RSIGS=RSIGS-SIGH(ICHN)  
        IF(RSIGS.LE.0.) GOTO 210    
  190   CONTINUE    
      ELSEIF(ISUB.EQ.96) THEN   
        CALL PYSPLI(MINT(11),21,KFL1,KFLDUM)    
        CALL PYSPLI(MINT(12),21,KFL2,KFLDUM)    
        MINT(1)=11  
        MINT(2)=1   
        IF(KFL1.EQ.KFL2.AND.RLU(0).LT.0.5) MINT(2)=2    
      ELSE  
        KFL1=21 
        KFL2=21 
        RSIGS=6.*RLU(0) 
        MINT(2)=1   
        IF(RSIGS.GT.1.) MINT(2)=2   
        IF(RSIGS.GT.2.) MINT(2)=3   
      ENDIF 
  210 IF(MINT(2).GT.10) THEN    
        MINT(1)=MINT(2)/10  
        MINT(2)=MOD(MINT(2),10) 
      ENDIF 
      MINT(15)=KFL1 
      MINT(16)=KFL2 
      MINT(13)=MINT(15) 
      MINT(14)=MINT(16) 
      VINT(141)=VINT(41)    
      VINT(142)=VINT(42)    
 1000 FORMAT(1X,'Error: maximum violated by',1P,E11.3,1X,   
     &'in event',1X,I7,'.'/1X,'Execution stopped!') 
 1100 FORMAT(1X,'ISUB = ',I3,'; Point of violation:'/1X,'tau=',1P, 
     &E11.3,', y* =',E11.3,', cthe = ',0P,F11.7,', tau'' =',1P,E11.3)   
      RETURN    
      END   
