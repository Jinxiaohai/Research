#!/bin/csh

    source /home/jhchen/.conenv.csh

    echo "====================================================="
    echo " Mission Started   at `date`"
    echo "====================================================="
    if ( $# != 1 ) then
    echo "====================================================="
    echo " Need One Parameter "
    echo " Program Exit      at `date`"
    echo "====================================================="
    exit
    endif
    #-------------------------------------------------
    set Folder = $PWD
    @ NBeg = 1
    @ NEnd = 1
    #-------------------------------------------------
    @ NCur=$NBeg
    while ( $NCur <= $NEnd )
	if ( $NCur > 0 && $NCur < 10 ) then
 		set PathId = /home/jinxiaohai/jinxiaohai/fluctuation/fluc/backup/ampt_0_$1
	# else if ( $NCur >= 10 && $NCur < 100 ) then
	# 	set PathId = $Folder/data/ampt_${1}_0$NCur
	# else if ( $NCur >= 100 && $NCur < 1000 ) then
        # set PathId = $Folder/data/ampt_${1}_$NCur
	else 
		echo "====================================================="
		echo " Exception ! NFile must < 1000 !"
		echo " Program Exit      at `date`"
		echo "====================================================="
		exit
	endif



        # generate data
    #_______________________________________________________________
        mkdir -p $PathId
        cp -r $Folder/ampt $PathId/
        cd $PathId/ampt
        make

        sh exec $1
    #_______________________________________________________________

        sleep 30


        create before fluc tree
    #_______________________________________________________________
        cp -r $Folder/zpcBeforeCreateTree $PathId
        find $PathId/ampt/ana/ -name zpc.dat > $PathId/zpcBeforeCreateTree/list/data.list
        cd $PathId/zpcBeforeCreateTree
        make clean 
        make
        cd $PathId
        set EXEH=$PathId/zpcBeforeCreateTree/bin/analysis
        $EXEH $PathId/zpcBeforeCreateTree/list/data.list $PathId/zpc_BeforeFluc
    #_______________________________________________________________
        
        sleep 30


        # create after fluc tree
    #_______________________________________________________________
        cp -r $Folder/zpcAfterCreateTree $PathId
        find $PathId/ampt/ana/ -name after_parton.dat > $PathId/zpcAfterCreateTree/list/data.list
        cd $PathId/zpcAfterCreateTree
        make clean 
        make
        cd $PathId
        set EXEH=$PathId/zpcAfterCreateTree/bin/analysis
        $EXEH $PathId/zpcAfterCreateTree/list/data.list $PathId/zpc_AfterFluc
    #_______________________________________________________________

        sleep 30

        # create ampt tree
    #_______________________________________________________________
        cp -r $Folder/amptCreateTree $PathId
        find $PathId/ampt/ana/ -name ampt.dat > $PathId/amptCreateTree/list/data.list
        cd $PathId/amptCreateTree
        make clean 
        make
        cd $PathId
        set EXEH=$PathId/amptCreateTree/bin/analysis
        $EXEH $PathId/amptCreateTree/list/data.list $PathId/ampt_AfterArt
    #_______________________________________________________________

        # get eccentricity
      #_______________________________________________________________
        cp -r $Folder/eccentricity    $PathId
        find $PathId  -name  zpc_BeforeFluc.root > $PathId/eccentricity/list/data.list
        cd $PathId/eccentricity
        make clean
        make
        cd $PathId
        set EXEH=$PathId/eccentricity/eccentricity.exe
        $EXEH  $PathId/eccentricity/list/data.list $PathId/eccentricity_BeforeFLuc.txt
      #_______________________________________________________________

        # sleep 30
      #_______________________________________________________________
        find $PathId  -name  zpc_AfterFluc.root > $PathId/eccentricity/list/data.list
        cd $PathId/eccentricity
        make clean
        make
        cd $PathId
        set EXEH=$PathId/eccentricity/eccentricity.exe
        $EXEH  $PathId/eccentricity/list/data.list $PathId/eccentricity_AfterFLuc.txt
        rm -rf $PathId/eccentricity
        rm -f saveFile.root
      #_______________________________________________________________

        # sleep 30
      #_______________________________________________________________
        cp -r $Folder/Eccentricity    $PathId
        find $PathId/ampt/ana -name npart-xy.dat > $PathId/Eccentricity/list/data.list
        cd $PathId/Eccentricity
        make clean
        make
        cd $PathId
        set EXEH=$PathId/Eccentricity/Eccentricity.exe
        $EXEH  $PathId/Eccentricity/list/data.list $PathId/Eccentricity_Initial.txt
        rm -rf $PathId/Eccentricity
      #_______________________________________________________________


	@ NCur += 1
    end

    echo "====================================================="
    echo " Mission Completed at `date`"
    echo "====================================================="
