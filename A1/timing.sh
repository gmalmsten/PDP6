echo \n >> timings.txt
echo New run >> timings.txt
for i in 1 2 4 8 16 32
do
  echo running on $i threads
  mpirun --bind-to none -np $i ./stencil ../../../maya/public/PDP_Assignment1/input100000000.txt out.txt 1 >> timings.txt
done
