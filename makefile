VARIANT := indexStruct
EXE := program
EXETYPE := opt
NUMPROCS := 1

main: setup cpuKernel gpuKernel
	mpicxx -O3 -I${PTL}/include -I./src -I/usr/local/cuda/include -DOPTL=3 -c src/main.cc -o obj/main.opt.o
	mpicxx -O0 -I${PTL}/include -I./src -I/usr/local/cuda/include -DOPTL=0 -c src/main.cc -o obj/main.dbg.o
	mpicxx -O3 ./obj/*.opt.o -o ${EXE}.opt -L${PTL}/lib -lPTL -L/usr/local/cuda/lib64 -lcudadevrt -lcudart
	mpicxx -O0 ./obj/*.dbg.o -o ${EXE}.dbg -L${PTL}/lib -lPTL -L/usr/local/cuda/lib64 -lcudadevrt -lcudart

run: main
	mpirun -np ${NUMPROCS} ./${EXE}.${EXETYPE} debug.ptl

cpuKernel:
	mpicxx -O3 -I${PTL}/include -I./src -I/usr/local/cuda/include -c src/Glob.cpp -o obj/Glob.opt.o
	mpicxx -O0 -I${PTL}/include -I./src -I/usr/local/cuda/include -c src/Glob.cpp -o obj/Glob.dbg.o
	mpicxx -O3 -I${PTL}/include -I./src -I/usr/local/cuda/include -c src/cpuKernel.cpp -o obj/cpuKernel.opt.o
	mpicxx -O0 -I${PTL}/include -I./src -I/usr/local/cuda/include -c src/cpuKernel.cpp -o obj/cpuKernel.dbg.o

gpuKernel:
	nvcc -ccbin=mpicxx -O3 -x cu -rdc=true -dc -I${PTL}/include -I./src -I/usr/local/cuda/include src/gpuKernel.cu -o obj/gpuKernel.opt.o
	nvcc -ccbin=mpicxx -O0 -x cu -rdc=true -dc -I${PTL}/include -I./src -I/usr/local/cuda/include src/gpuKernel.cu -o obj/gpuKernel.dbg.o
	nvcc -ccbin=mpicxx -rdc=true -dlink obj/gpuKernel.opt.o -o obj/gpuKernel.dlink.opt.o
	nvcc -ccbin=mpicxx -rdc=true -dlink obj/gpuKernel.dbg.o -o obj/gpuKernel.dlink.dbg.o

setup:
	unlink src || echo "no src"
	ln -sf variations/${VARIANT} src
	mkdir -p obj
	mkdir -p output

clean:
	unlink src || echo "no src"
	rm -f ${EXE}.*
	rm -rf obj
	rm -rf output
