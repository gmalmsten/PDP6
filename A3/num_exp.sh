rm weak.txt
echo Weak scaling >> weak.txt
echo 1 >> weak.txt
mpirun -np 1 qsort 125000000 out.txt 1 >> weak.txt
mpirun -np 2 qsort 250000000 out.txt 1 >> weak.txt
mpirun -np 4 qsort 500000000 out.txt 1 >> weak.txt
mpirun -np 8 qsort 1000000000 out.txt 1 >> weak.txt
mpirun -np 16 qsort 2000000000 out.txt 1 >> weak.txt

echo 2 >> weak.txt
mpirun -np 1 qsort 125000000 out.txt 2 >> weak.txt
mpirun -np 2 qsort 250000000 out.txt 2 >> weak.txt
mpirun -np 4 qsort 500000000 out.txt 2 >> weak.txt
mpirun -np 8 qsort 1000000000 out.txt 2 >> weak.txt
mpirun -np 16 qsort 2000000000 out.txt 2 >> weak.txt

echo 3 >> weak.txt
mpirun -np 1 qsort 125000000 out.txt 3 >> weak.txt
mpirun -np 2 qsort 250000000 out.txt 3 >> weak.txt
mpirun -np 4 qsort 500000000 out.txt 3 >> weak.txt
mpirun -np 8 qsort 1000000000 out.txt 3 >> weak.txt
mpirun -np 16 qsort 2000000000 out.txt 3 >> weak.txt


