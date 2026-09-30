mpic++ -c mpmd.cpp -o mpmd.o -fPIC -O3 -fopenmp
ar rcs ../lib/libmpmd.a mpmd.o

# mpic++ proa.cpp -o proa -L./ -lmpmd -O3 -fopenmp

# mpic++ prob.cpp -o prob -L./ -lmpmd -O3 -fopenmp
# mpirun -n 2 ./proa : -n 2 ./prob