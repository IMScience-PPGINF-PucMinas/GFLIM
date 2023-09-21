BIN=./bin
DEMO=./demo
INCLUDE=./include
OBJ=./obj
SRC=./src


FLAGS = -fPIC -O3 -Wall -lm
INCLUDES = -I $(INCLUDE)


all: $(OBJ)/znzlib.o $(OBJ)/nifti1_io.o

demos: clib_01_read_write


# generic src compilation
$(OBJ)/%.o: $(SRC)/%.c $(INCLUDE)/%.h
	$(CC) $(FLAGS) -c $< $(INCLUDES) -o $@ 


%: $(DEMO)/%.c
	$(CC) $< $(OBJ)/*.o $(INCLUDES) -o $(BIN)/$@ $(FLAGS)



clean:
	$(RM) $(OBJ)/*

