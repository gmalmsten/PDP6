CC = mpicc
CCFLAGS = -g -O3
LIBS = -lm

BINS = q_sort

all: $(BINS)

%: %.c
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(BINS)

