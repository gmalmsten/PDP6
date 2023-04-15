###############################################################################
# Makefile for assignment 1, Parallel and Distributed Computing 2022.
###############################################################################

CC = mpicc
CFLAGS = -std=c99 -g -Ofast -Wall 
LIBS = -lm

BIN = matmul

all: $(BIN)

matmul: matmul.c 
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)


	
clean:
	$(RM) $(BIN)
