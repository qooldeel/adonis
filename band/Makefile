EXEC = test

MAIN_PROG = main.o

OBJS = $(MAIN_PROG)   #plus possible further objects

HEAD_DIR =./

ERR = -Wfatal-errors

CC_FLAGS = -Wall -O3 $(ERR) $(DEFS) -pedantic

CC = c++

INCL = -I$(HEAD_DIR) 


.cc.o:
	@echo ">>>> COMPILING..."
	$(CC) $(CC_FLAGS) $(INCL)  -c $*.cc

$(EXEC): $(OBJS) 
	@echo ">>>> LINKING..."
	$(CC) $(CC_FLAGS) -o $(EXEC) $(OBJS)
	@echo "Done! :) :) :)"

#clean up at the very end
.PHONY: clean
clean:
	rm -rf $(OBJS) $(EXEC) *.out
