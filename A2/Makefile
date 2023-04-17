###############################################################################
# Makefile for assignment 2, Parallel and Distributed Computing 2022.
###############################################################################

CC = mpicc
CFLAGS = -std=c99 -g -Ofast -Wall 
LIBS = -lm -lblas

BIN = matmul

all: $(BIN)

matmul: matmul.c 
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)


clean:
	$(RM) $(BIN)
