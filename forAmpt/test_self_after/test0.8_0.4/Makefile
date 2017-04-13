##########################################################
##                 Created  by  xiaohai                 ##
##                 Telphone : 18501781924               ##
##            E-mail : jinxiaohai@sinap.ac.cn           ##
##            E-mail : xiaohaijin@outlook.com           ##
##   Address : Shanghai Institute of Applied Physics    ##
##########################################################

EXE = zpc.exe
OBJ = zpc.o
ROOTLIBS     =   $(shell root-config --libs)
#ROOTLIBS  =  -L/home/xiaohai/SoftWear/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript-lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic
ROOTGLIBS    =   $(shell root-config --glibs)
#ROOTGLIBS  =  -L/home/xiaohai/SoftWear/root/lib -lGui -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic

INCDIR = -I$(ROOTSYS)/include
LIBDIR = -L$(ROOTSYS)/lib
CC = g++
FLAGS = -Wall -g $(INCDIR) $(LIBDIR)
COMPILE = $(CC) $(FLAGS) -c

$(EXE):$(OBJ)
	$(CC) -o $(EXE) $(OBJ) $(ROOTLIBS) $(ROOTGLIBS)

$(OBJ):zpc.cpp head.h   ;       $(COMPILE) zpc.cpp

.PHONY:clean
clean:
	-rm $(OBJ) $(EXE)
