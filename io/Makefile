EXEC = test

OBJS = main.o 

INCL = -I$./

CC = c++

WFATAL = -Wfatal-errors

OPTIMIZE = -O3 #-funroll-loops

CC_FLAGS = -Wall -g $(OPTIMIZE) $(WFATAL)

.cc.o:
	@echo "Starting compilation, Marc..."
	$(CC) $(CC_FLAGS) $(INCL) -c $*.cc 

$(EXEC): $(OBJS)
	@echo "Linking and building executable..."
	$(CC) -o $(EXEC) $(OBJS)

.PHONY: clean
clean:
	rm -rf $(OBJS) $(EXEC) *.out
