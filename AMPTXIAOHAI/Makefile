# was created by xiaohai.

F77 = gfortran

## source file path
SRC_PATH   := ./src
SRCS         += $(wildcard $(SRC_PATH)/*.f)
OBJS        := $(SRCS:.f=.o)

all: clean build

.PHONY.:build
build:
	@$(F77) $(SRCS) -o ampt

.PHONY.:clean
clean:
	@$(RM) $(OBJS) ampt  
