rm backwards.txt
echo Weak scaling >> backwards.txt
echo 1 >> backwards.txt
mpirun -np 1 qsort_b 125000000 out.txt 1 >> backwards.txt
mpirun -np 2 qsort_b 250000000 out.txt 1 >> backwards.txt
mpirun -np 4 qsort_b 500000000 out.txt 1 >> backwards.txt
mpirun -np 8 qsort_b 1000000000 out.txt 1 >> backwards.txt
mpirun -np 16 qsort_b 2000000000 out.txt 1 >> backwards.txt

echo 2 >> backwards.txt
mpirun -np 1 qsort_b 125000000 out.txt 2 >> backwards.txt
mpirun -np 2 qsort_b 250000000 out.txt 2 >> backwards.txt
mpirun -np 4 qsort_b 500000000 out.txt 2 >> backwards.txt
mpirun -np 8 qsort_b 1000000000 out.txt 2 >> backwards.txt
mpirun -np 16 qsort_b 2000000000 out.txt 2 >> backwards.txt

echo 3 >> backwards.txt
mpirun -np 1 qsort_b 125000000 out.txt 3 >> backwards.txt
mpirun -np 2 qsort_b 250000000 out.txt 3 >> backwards.txt
mpirun -np 4 qsort_b 500000000 out.txt 3 >> backwards.txt
mpirun -np 8 qsort_b 1000000000 out.txt 3 >> backwards.txt
mpirun -np 16 qsort_b 2000000000 out.txt 3 >> backwards.txt


echo Strong scaling >> backwards.txt
echo 1 >> backwards.txt
mpirun -np 1 qsort_b 250000000 out.txt 1 >> backwards.txt
mpirun -np 2 qsort_b 250000000 out.txt 1 >> backwards.txt
mpirun -np 4 qsort_b 250000000 out.txt 1 >> backwards.txt
mpirun -np 8 qsort_b 250000000 out.txt 1 >> backwards.txt
mpirun -np 16 qsort_b 250000000 out.txt 1 >> backwards.txt

echo 2 >> backwards.txt
mpirun -np 1 qsort_b 250000000 out.txt 2 >> backwards.txt
mpirun -np 2 qsort_b 250000000 out.txt 2 >> backwards.txt
mpirun -np 4 qsort_b 250000000 out.txt 2 >> backwards.txt
mpirun -np 8 qsort_b 250000000 out.txt 2 >> backwards.txt
mpirun -np 16 qsort_b 250000000 out.txt 2 >> backwards.txt

echo 3 >> backwards.txt
mpirun -np 1 qsort_b 250000000 out.txt 3 >> backwards.txt
mpirun -np 2 qsort_b 250000000 out.txt 3 >> backwards.txt
mpirun -np 4 qsort_b 250000000 out.txt 3 >> backwards.txt
mpirun -np 8 qsort_b 250000000 out.txt 3 >> backwards.txt
mpirun -np 16 qsort_b 250000000 out.txt 3 >> backwards.txt
