SRC = ./
BIN = ./

# CC should be set to the name of your favorite C compiler.
CC = cc

# Where to find the basic X Windows (X11) #include files and libraries
SHOWMESWITCHES = -I/usr/X11R6/include -L/usr/X11R6/lib

default: $(SRC)visualization.c
	$(CC) $(CSWITCHESFAST) $(SHOWMESWITCHES) -o $(BIN)visualization \
		$(SRC)visualization.c -lX11 -lm


clean:
	-$(RM) $(BIN)visualization *~
