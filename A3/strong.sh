echo Strong scaling >> strong.txt
echo 1 >> strong.txt
mpirun -np 1 qsort 250000000 out.txt 1 >> strong.txt
mpirun -np 2 qsort 250000000 out.txt 1 >> strong.txt
mpirun -np 4 qsort 250000000 out.txt 1 >> strong.txt
mpirun -np 8 qsort 250000000 out.txt 1 >> strong.txt
mpirun -np 16 qsort 250000000 out.txt 1 >> strong.txt

echo 2 >> strong.txt
mpirun -np 1 qsort 250000000 out.txt 2 >> strong.txt
mpirun -np 2 qsort 250000000 out.txt 2 >> strong.txt
mpirun -np 4 qsort 250000000 out.txt 2 >> strong.txt
mpirun -np 8 qsort 250000000 out.txt 2 >> strong.txt
mpirun -np 16 qsort 250000000 out.txt 2 >> strong.txt

echo 3 >> strong.txt
mpirun -np 1 qsort 250000000 out.txt 3 >> strong.txt
mpirun -np 2 qsort 250000000 out.txt 3 >> strong.txt
mpirun -np 4 qsort 250000000 out.txt 3 >> strong.txt
mpirun -np 8 qsort 250000000 out.txt 3 >> strong.txt
mpirun -np 16 qsort 250000000 out.txt 3 >> strong.txt