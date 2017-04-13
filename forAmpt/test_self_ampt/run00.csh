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
 		# set PathId = /home/jinxiaohai/xiaohai/test_zpc/
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
	# mkdir -p $PathId

	# cp -r $Folder/test0.8_0.4/ $PathId/
#    cp -r $Folder/createTree $PathId/
	cd $Folder/test0.8_0.4
        # find /home/jinxiaohai/scratch/fluctuation/fluc/b0 -name ampt.dat  > ./list/data.list
	# make
        ./run.sh
# 	@ nrandom = `date '+%s%N'` + $1
# #    @ nrandom = `date '+%d%H%M%S'` + $1
# 	sh exec $nrandom
    
# 	mv $PathId/ampt/ana $PathId/
#         mv $PathId/ampt/hiji*.f $PathId/
# 	rm -rf $PathId/ampt
#	cd $PathId/createTree
#make clean 
#	make
#    cd $PathId
#    set EXEH=$PathId/createTree/bin/analysis
#    $EXEH $PathId/ampt/ana/ampt.dat $PathId/ampt_afterART	
	
	@ NCur += 1
end

echo "====================================================="
echo " Mission Completed at `date`"
echo "====================================================="
