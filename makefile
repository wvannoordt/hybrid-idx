EXE := program
NUMPROCS := 6

main: setup cpuKernel gpuKernel
	mpicxx -I${PTL}/include -I./src -c src/main.cc -o obj/main.o
	mpicxx ./obj/*.o -o ${EXE} -L${PTL}/lib -lPTL

run: main
	mpirun -np ${NUMPROCS} ./${EXE}

cpuKernel:
	@echo "(no cpu yet)"

gpuKernel:
	@echo "(no gpu yet)"

setup:
	mkdir -p obj

clean:
	rm -f ${EXE}
	rm -rf obj
