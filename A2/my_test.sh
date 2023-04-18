#!/bin/bash
np=4

make
rm -f out.txt
mpirun -np 2 ./matmul input${np}.txt out.txt
echo "---Printing diff---"
diff -s out.txt ref_input${np}.txt