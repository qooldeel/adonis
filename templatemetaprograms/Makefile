#================= Change me, Marc ==========================================

EXEC = testmeta

MAIN_PROG = testmeta.o

OBJS = $(MAIN_PROG) 

HEAD_DIR =./  #include subdirectory stuff (if any)
CPPAD = ${HOME}/SOFTWARE/AutoDiff
INCL = -I$(HEAD_DIR) -I$(CPPAD)

LIBS = -llapack -lblas

WFAT =  -Wfatal-errors

DEFS = -DLAPACK_USAGE


#==========================================================================

CC = c++

CC_FLAGS = -O3 -Wall -pg -funroll-loops -finline-functions -DNDEBUG $(WFAT) $(DEFS)


#specify compilation of C++ source files to object files (u can use $< instead of $*.cc in the following block.  
.cc.o:
	@echo "Compling..."
	$(CC) $(CC_FLAGS) $(INCL) -c $*.cc 



$(EXEC): $(OBJS)
	@echo "Writing executable \"$(EXEC)\", Marc..."
	$(CC) -o $(EXEC) $(OBJS) $(LIBS) 


#clean up at the very end
.PHONY: clean
clean:
	rm -rf $(OBJS) $(EXEC) *.out
