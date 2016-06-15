
###  USAGE:
## 1.) replace FILENAME = with the name of the file you wanna to compile 
## 2.) go back to terminal and execute: make clean; make; ./file
##     where 'file' is the executable from EXEC 

## ========== TODO: replace FILENAME (without suffix) by specific name of file
##                  you're going to compile, e.g. vec, blatt6, etc.
##  CAUTION: BEWARE OF WHITE CHARS IN FILENAME!!! (otherwise your file might
##           not be found ;) ) 

FILENAME = odetest
#FILENAME = odetest

## Quickstart
## Instructions on how to quickly get started with lcov:

## Ensure that the project is built using GCC
## Add --coverage to compiler and linker flags (for example CFLAGS and LDFLAGS)
## Compile and run
## Collect coverage data:
## lcov --capture --directory project-dir --output-file coverage.info
## Generate HTML output:
## genhtml coverage.info --output-directory out

## enable code coverage (uncomment RHS for no coverage)
COVERAGE = #--coverage  
##============================================================================



EXEC = $(FILENAME)

MAIN_PROG = $(FILENAME).o

OBJS = $(MAIN_PROG)   #plus possible further objects

## uncomment the line starting with 'PROFILE'. Then compile your prog. and
## execute it (possibly with arguments) as usual. Profile data is written
## to `gmon.out' just before exiting (all files with the same name are going
## to be overwritten). Profile data from `gmon.out' is then interpreted by 
## executing e.g.
##            gprof ./executable > myprofile.txt
## By editing `myprofile.txt', functions, members, etc. 
## with decreasing execution time are listed.
## NOTE: profiling can last a considerable time! So don't get nervous... 
#PROFILE = -pg -O3 -funroll-loops

#if you have some/all header files in a subdirectory 
HEAD_DIR =./

THRDBASE = $(HOME)/SOFTWARE

#CPPAD = $(THRDBASE)/AutoDiff

### BEGIN:  UMFPACK-SPECIFIC STUFF ####
UMFPACK_INCL = -I$(THRDBASE)/UMFPACK/Include -I$(THRDBASE)/AMD/Include -I$(THRDBASE)/SuiteSparse_config

UMFPACK_CONFIG = -DNCHOLMOD #-DLP64

UMFPACK_LIBS = $(THRDBASE)/UMFPACK/Lib/libumfpack.a $(THRDBASE)/AMD/Lib/libamd.a $(THRDBASE)/SuiteSparse_config/libsuitesparseconfig.a
### END: UMFPACK-SPECIFIC STUFF ####



INCL = -I$(HEAD_DIR) $(UMFPACK_INCL) #-I$(CPPAD)
LIB = -lm -lrt
# external libraries, for instance BLAS, LAPACK, ATLAS 
LIBS = -llapack -lblas -lgfortran $(LIB) 

LIBS_ADD = 

# linking libraries
LCC_LIBS = -lm #-lg2c

# flags for linking
LCC_FLAGS = #-L/usr/lib/liblapack.so.3gf -L/usr/lib/sse2

LCC_FLAGS_ADD = 

BUG = #-DNDEBUG

DEFS = -DLAPACK_USAGE 

#my own libraries (if any exist so far)
MY_LIBS =

PEDANTIC = -pedantic

#show only first occuring error
ERR = -Wfatal-errors $(PEDANTIC) 

CC = c++


#Uncomment the '-D...' if you don't want to compile the stuff
ADD_CC_FLAGS = -DDUNE_STUFF   

WALL = -Wall

#CC_FLAGS = -O3 -Wall -g -pg $(ADD_CC_FLAGS) $(DEFS) $(BUG) $(ERR)
CC_FLAGS =  $(WALL) $(ADD_CC_FLAGS) $(DEFS) $(BUG) $(ERR) $(COVERAGE)


#specify compilation of C++ source files to object files (u can use $< instead of $*.cc in the following block.  
.cc.o:
	@echo "Starting compilation, Marc..."
	$(CC) $(CC_FLAGS) $(PROFILE) $(UMFPACK_CONFIG) $(INCL) -c $*.cc 



$(EXEC): $(OBJS)
	@echo "Linking and building executable..."
	$(CC) $(CC_FLAGS) $(UMFPACK_CONFIG) -o $(EXEC) $(OBJS) $(UMFPACK_LIBS) $(PROFILE) $(LCC_FLAGS) $(LCC_FLAGS_ADD) $(LIBS) $(LIBS_ADD) $(LCC_LIBS) 
	@echo "Done!"

#clean up at the very end
.PHONY: clean
clean:
	rm -rf $(OBJS) $(EXEC) *.out
