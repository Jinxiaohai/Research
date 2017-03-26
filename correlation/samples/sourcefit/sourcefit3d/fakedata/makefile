CORALHOMEDIR = ../../../../../coralpha
CPP = g++ -bind_at_load
OPT = ${CPPFLAGS}

# You shouldn't have to adjust anything below
THISDIR = ${PWD}

coral :
	cd ${CORALHOMEDIR}/lib; make; cd ${THISDIR}

fakedatawriter_3dgaussian : fakedatawriter_3dgaussian.cc
	${CPP} ${OPT} -L${CORALHOMEDIR}/lib -lgsl -lgslcblas -lcoral -I ${CORALHOMEDIR}/include fakedatawriter_3dgaussian.cc -o fakedatawriter_3dgaussian

fakedatawriter_blast : fakedatawriter_blast.cc 
	${CPP} ${OPT} -L${CORALHOMEDIR}/lib -lgsl -lgslcblas -lcoral -I ${CORALHOMEDIR}/include fakedatawriter_blast.cc -o fakedatawriter_blast
