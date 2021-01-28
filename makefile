EXE := program
NUMPROCS := 1

main: setup cpuKernel gpuKernel
	mpicxx -I${PTL}/include -I./src -I/usr/local/cuda/include -c src/main.cc -o obj/main.o
	mpicxx ./obj/*.o -o ${EXE} -L${PTL}/lib -lPTL -L/usr/local/cuda/lib64 -lcudadevrt -lcudart

run: main
	mpirun -np ${NUMPROCS} ./${EXE}

cpuKernel:
	mpicxx -I${PTL}/include -I./src -I/usr/local/cuda/include -c src/Glob.cpp -o obj/Glob.o
	mpicxx -I${PTL}/include -I./src -I/usr/local/cuda/include -c src/cpuKernel.cpp -o obj/cpuKernel.o

gpuKernel:
	nvcc -ccbin=mpicxx -O0 -x cu -rdc=true -dc -I${PTL}/include -I./src -I/usr/local/cuda/include src/gpuKernel.cu -o obj/gpuKernel.o
	nvcc -ccbin=mpicxx -rdc=true -dlink obj/gpuKernel.o -o obj/gpuKernel.dlink.o

setup:
	mkdir -p obj

clean:
	rm -f ${EXE}
	rm -rf obj
