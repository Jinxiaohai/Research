CORALHOME=..
CPP = /usr/bin/g++ -bind_at_load
OPT=${CFLAGS}
#
CORALINCLUDE=${CORALHOME}/include
CORALLIB=${CORALHOME}/lib
THISIDIR=${PWD}
coral :
	cd ${CORALLIB}; make; cd ${THISDIR}
wfsample : wfsample.cc ${CORALLIB}/libcoral.a
	${CPP} ${OPT} -L${CORALLIB} -lgsl -lgslcblas -lcoral -I ${CORALINCLUDE} wfsample.cc -o wfsample
kernelsample : kernelsample.cc  ${CORALLIB}/libcoral.a
	${CPP} ${OPT} -L${CORALLIB} -lgsl -lgslcblas -lcoral -I ${CORALINCLUDE} kernelsample.cc -o kernelsample
sourcesample_gauss : sourcesample_gauss.cc ${CORALLIB}/libcoral.a
	${CPP} ${OPT} -L${CORALLIB} -lgsl -lgslcblas -lcoral -I ${CORALINCLUDE} sourcesample_gauss.cc -o sourcesample_gauss
sourcesample_OSCAR : sourcesample_OSCAR.cc  ${CORALLIB}/libcoral.a
	${CPP} ${OPT} -L${CORALLIB} -lgsl -lgslcblas -lcoral -I ${CORALINCLUDE} sourcesample_OSCAR.cc -o sourcesample_OSCAR
sourcesample_blast : sourcesample_blast.cc ${CORALLIB}/libcoral.a
	${CPP} ${OPT} -L${CORALLIB} -lgsl -lgslcblas -lcoral -I ${CORALINCLUDE} sourcesample_blast.cc -o sourcesample_blast
lisample : lisample.cc ${CORALLIB}/libcoral.a
	${CPP} ${OPT} -L${CORALLIB} -lgsl -lgslcblas -lcoral -I ${CORALINCLUDE} lisample.cc -o lisample
