#include "gpuKernel.h"
#include <iostream>
__global__ void K_TestFunctionsGpu(double* flow, const InputClass& input)
{
    
}

void TestFunctionsGpu(double* flow, const InputClass& input)
{
    //Global memory:            7.907288 GB
    //Shared memory per block:  48.000000 KB
    //Warp size:                32￼
    //Max threads per block:    1024
    //Max thread dimension:     1,024  x  1,024  x  64
    //Max grid size:            2,147,483,647  x  65,535  x  65,535
    //Total constant memory:    64.000000 KB
    
    dim3 blockConf;
    blockConf.x = BLOCK_SIZE;
    blockConf.y = BLOCK_SIZE;
#if(IS3D)
    blockConf.z = BLOCK_SIZE;
#endif
    dim3 gridConf;
    int numcells[DIM];
    for (int i = 0; i < DIM; i++) {numcells[i] = input.nxb[i] + 2*input.nguard;}
    gridConf.x = input.lnblocks*((numcells[0] + BLOCK_SIZE - 1)/BLOCK_SIZE);
    gridConf.y = (numcells[1] + BLOCK_SIZE - 1)/BLOCK_SIZE;
#if(IS3D)
    gridConf.z = (numcells[2] + BLOCK_SIZE - 1)/BLOCK_SIZE;
#endif

    if (mypenoG==0)
    {
        std::cout << "GP Config:\nblock: " << blockConf.x << " x " << blockConf.y;
        if (IS3D) std::cout << " x " << blockConf.z;
        std::cout << "\ngrid:  " << gridConf.x << " x " << gridConf.y;
        if (IS3D) std::cout << " x " << gridConf.z;
        std::cout << std::endl;
    }
}