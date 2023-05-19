CC = mpicc
CCFLAGS = -g -Ofast -march=native
LIBS = -lm

BINS = qsort

all: $(BINS)

%: %.c
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(BINS)

