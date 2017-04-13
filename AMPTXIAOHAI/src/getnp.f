        subroutine getnp
        PARAMETER (MAXSTR=150001)
        COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
        COMMON/HMAIN2/KATT(MAXSTR,4),PATT(MAXSTR,4)
        COMMON /HPARNT/HIPR1(100), IHPR2(50), HINT1(100), IHNT2(50)
        common/snn/efrm,npart1,npart2,epsiPz,epsiPt,PZPROJ,PZTARG
        SAVE   
        if(NATT.eq.0) then
           npart1=0
           npart2=0
           return
        endif
        PZPROJ=SQRT(HINT1(6)**2-HINT1(8)**2)
        PZTARG=SQRT(HINT1(7)**2-HINT1(9)**2)
        epsiPz=0.01
        epsiPt=1e-6
        nspec1=0
        nspec2=0
        DO 1000 I = 1, NATT
           if((KATT(I,1).eq.2112.or.KATT(I,1).eq.2212)
     1          .and.abs(PATT(I, 1)).le.epsiPt
     2          .and.abs(PATT(I, 2)).le.epsiPt) then
              if(PATT(I, 3).gt.amax1(0.,PZPROJ-epsiPz)) then
                 nspec1=nspec1+1
              elseif(PATT(I, 3).lt.(-PZTARG+epsiPz)) then
                 nspec2=nspec2+1
              endif
           endif
 1000   CONTINUE
        npart1=IHNT2(1)-nspec1
        npart2=IHNT2(3)-nspec2
        return
        end
