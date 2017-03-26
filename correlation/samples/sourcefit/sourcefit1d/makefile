CMINUITHOMEDIR = ../../../../cminuit
CORALHOMEDIR = ../../../../coralpha
CPP = g++ -bind_at_load
OPT = ${CPPFLAGS}
MINUITFLAGS = -bind_at_load -Df2cFortran -lfrtbegin -lg2c -lSystemStubs -lSystem -lm

# You shouldn't have to adjust anything below
THISDIR = ${PWD}
sourcefit_minuit_GX1d : sourcefit_minuit_GX1d.cc ${CORALHOMEDIR}/lib/libcoral.a
	${CPP} ${MINUITDEFS} ${OPT} ${MINUITFLAGS} sourcefit_minuit_GX1d.cc -I ${CMINUITHOMEDIR}/include -I ${CORALHOMEDIR}/include ${CMINUITHOMEDIR}/lib/libcminuit.a ${CORALHOMEDIR}/lib/libcoral.a -lgsl -lgslcblas -o sourcefit_minuit_GX1d

sourcefit_GX1d : sourcefit_GX1d.cc ${CORALHOMEDIR}/lib/libcoral.a
	${CPP} ${MINUITDEFS} ${OPT} sourcefit_GX1d.cc -I ${CMINUITHOMEDIR}/include -I ${CORALHOMEDIR}/include ${CORALHOMEDIR}/lib/libcoral.a -lgsl -lgslcblas -o sourcefit_GX1d

coral :
	cd ${CORALHOMEDIR}/lib; make; cd ${THISDIR}

cminuit :
	cd ${CMINUITHOMEDIR}/lib; make libcminuit.a; cd ${THISDIR};

