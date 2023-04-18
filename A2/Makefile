###############################################################################
# Makefile for assignment 2, Parallel and Distributed Computing 2022.
###############################################################################

CC = mpicc
CFLAGS = -std=c99 -g -Ofast 
LIBS = -lm -lopenblas

BIN = matmul matmul_backup

all: $(BIN)

matmul: matmul.c 
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

matmul_backup: matmul_backup.c 
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(BIN)
