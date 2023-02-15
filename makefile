CC = gcc
CFLAGS = -lpthread -lm -O3
OBJ= 1_fetch.o 2_queue.o 3_inspect.o 4_ode.o 5_equation.o

all : trnode
trnode: $(OBJ)
	$(CC) $(OBJ) -o $@ $(CFLAGS)
	rm *.o
$(OBJ): %.o: %.c 0_header.h
	$(CC) -c $< $(CFLAGS)

clear:
	rm -rf trnode debug *.o data/20*
debug:
	gcc 0_header.h 1_fetch.c 2_queue.c 3_inspect.c 4_ode.c 5_equation.c -o debug -lpthread -lm -g -O0