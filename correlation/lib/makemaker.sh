#!/bin/sh
export CPP="g++ -bind_at_load"
echo CPP=${CPP}
echo OPT=\$\{CFLAGS\}
echo \#
echo \# You shouldn\'t have to mess with much below here
export CORALSRC="../src"
export CORALINC="../include"
export HCFILES=""

echo CORALSRC=${CORALSRC}
echo CORALINC=${CORALINC}
for ifile in `find ${CORALSRC} \( -name \*.h -o -name \*.cc \)`
do
    #echo ifile=${ifile}
    length=`echo ${ifile} | awk '{print length ($1)}'`
    n=`expr ${length} - 3`
    last=0
    while [ ${last} = 0 ]
    do
	char=`echo ${ifile} ${n} | awk '{print substr($1,$2,1)}'`
	#echo ${char}
	if [ ${char} = "/" ]
	then
	    last=${n}
	fi
	n=`expr ${n} - 1`
    done
    #echo length of ${ifile} = ${length}
    #echo last / is at ${last}
    jfile=`echo ${ifile} ${last} | awk '{print substr($1,$2+1,length($1)-$2)}'`
    HCFILES=`echo ${HCFILES} ${CORALINC}/${jfile}`
done
echo HCFILES=${HCFILES}
echo \#
echo libcoral.a : coral.o
printf "\tar -ru libcoral.a coral.o\n"
echo coral.o : \$\{CORALINC\}/coral.cc \$\{HCFILES\}
printf "\t\${CPP} -c -FPIC -I\${CORALINC} \${OPT} -o coral.o \${CORALINC}/coral.cc\n"
echo clean :
printf "\trm -f *.o *~\n"
echo \#
for ifile in `find ${CORALSRC} \( -name \*.h -o -name \*.cc \)`
do
    #echo ifile=${ifile}
    length=`echo ${ifile} | awk '{print length ($1)}'`
    n=`expr ${length} - 3`
    last=0
    while [ ${last} = 0 ]
    do
	char=`echo ${ifile} ${n} | awk '{print substr($1,$2,1)}'`
	#echo ${char}
	if [ ${char} = "/" ]
	then
	    last=${n}
	fi
	n=`expr ${n} - 1`
    done
    #echo length of ${ifile} = ${length}
    #echo last / is at ${last}
    jfile=`echo ${ifile} ${last} | awk '{print substr($1,$2+1,length($1)-$2)}'`
    echo ${CORALINC}/${jfile} : ${ifile}
    printf "\tcp $ifile ${CORALINC}/\n"
done
