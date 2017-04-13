        FUNCTION FNSTRU(X)
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
        SAVE   
        FNSTRU=(1.0-X)**HIPR1(44)/
     &                (X**2+HIPR1(45)**2/HINT1(1)**2)**HIPR1(46)
        RETURN
        END
