CC = mpicc
CCFLAGS = -g -Ofast
LIBS = -lm

BINS = quicksort

all: $(BINS)

%: %.c
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(BINS)

