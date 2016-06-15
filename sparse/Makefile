
EXEC = main

MAIN_PROG = main.o

OBJS = $(MAIN_PROG)   #plus possible further objects

HEAD_DIR =./

ERR = -Wfatal-errors

DEFS = -DLAPACK_USAGE   # apply Lapack

CC_FLAGS = -Wall -O3 $(ERR) $(DEFS) -pedantic

CC = c++


### BEGIN:  UMFPACK-SPECIFIC STUFF ####
THRDBASE = $(HOME)/SOFTWARE

UMFPACK_INCL = -I$(THRDBASE)/UMFPACK/Include -I$(THRDBASE)/AMD/Include -I$(THRDBASE)/SuiteSparse_config

UMFPACK_CONFIG = -DNCHOLMOD #-DLP64

UMFPACK_LIBS = $(THRDBASE)/UMFPACK/Lib/libumfpack.a $(THRDBASE)/AMD/Lib/libamd.a $(THRDBASE)/SuiteSparse_config/libsuitesparseconfig.a
### END: UMFPACK-SPECIFIC STUFF ####


INCL = -I$(HEAD_DIR) $(UMFPACK_INCL)

LIB = -lm -lrt
LIBS = -lblas -llapack -lgfortran $(LIB) 


.cc.o:
	@echo ">>>> COMPILING..."
	$(CC) $(CC_FLAGS) $(UMFPACK_CONFIG)  $(INCL)  -c $*.cc

$(EXEC): $(OBJS) 
	@echo ">>>> LINKING..."
	$(CC) $(CC_FLAGS) $(UMFPACK_CONFIG) -o $(EXEC) $(OBJS) $(UMFPACK_LIBS) $(LIBS)
	@echo "Done! :) :) :)"

#clean up at the very end
.PHONY: clean
clean:
	rm -rf $(OBJS) $(EXEC) *.out
