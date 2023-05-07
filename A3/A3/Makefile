CC = mpicc
CCFLAGS = -g -Ofast -march=native
LIBS = -lm

BINS = quicksort

all: $(BINS)

%: %.c
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(BINS)

