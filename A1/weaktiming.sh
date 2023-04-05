for i in 1 2 4 8
do
  echo running on $i threads
  mpirun --bind-to none -np $i ./stencil ../../../maya/public/PDP_Assignment1/input${i}000000.txt out.txt 1 >> timetest.txt
done