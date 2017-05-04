##########################################################
##                 Created  by  xiaohai                 ##
##                 Telphone : 18501781924               ##
##            E-mail : jinxiaohai@sinap.ac.cn           ##
##            E-mail : xiaohaijin@outlook.com           ##
##   Address : Shanghai Institute of Applied Physics    ##
##########################################################

EXE = partoninitial.exe
OBJ = partoninitial.o
ROOTLIBS     =   $(shell root-config --libs)
#ROOTLIBS  =  -L/home/xiaohai/SoftWear/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript-lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic
ROOTGLIBS    =   $(shell root-config --glibs)
#ROOTGLIBS  =  -L/home/xiaohai/SoftWear/root/lib -lGui -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic

INCDIR = -I$(ROOTSYS)/include
LIBDIR = -L$(ROOTSYS)/lib
CC = g++
FLAGS = -Wall -g -std=c++11 $(INCDIR) $(LIBDIR)
COMPILE = $(CC) $(FLAGS) -c

$(EXE):$(OBJ)
	$(CC) -o $(EXE) $(OBJ) $(ROOTLIBS) $(ROOTGLIBS)

partoninitial.o:partoninitial.cpp  ;       $(COMPILE) partoninitial.cpp

.PHONY:clean
clean:
	-rm $(OBJ) $(EXE)
