CORALHOMEDIR = ../../../../../coralpha
CPP = g++ -bind_at_load
OPT = ${CPPFLAGS}

# You shouldn't have to adjust anything below
THISDIR = ${PWD}

coral :
	cd ${CORALHOMEDIR}/lib; make; cd ${THISDIR}

kwriter : kwriter.cc 
	${CPP} ${OPT} -L${CORALHOMEDIR}/lib -lgsl -lgslcblas -lcoral -I ${CORALHOMEDIR}/include kwriter.cc -o kwriter

fakedatawriter_GX1d : fakedatawriter_GX1d.cc
	${CPP} ${OPT} -L${CORALHOMEDIR}/lib -lgsl -lgslcblas -lcoral -I ${CORALHOMEDIR}/include fakedatawriter_GX1d.cc -o fakedatawriter_GX1d
