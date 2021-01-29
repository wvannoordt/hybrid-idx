#include "gpuKernel.h"
#include "Idx.h"
#include "CuErr.h"
#include <iostream>
__global__ void K_TestFunctionsGpu(double* flow, const InputClass input, const int lb)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x - input.nguard;
    int j = threadIdx.y + blockIdx.y*blockDim.y - input.nguard;
#if(IS3D)    
    int k = threadIdx.z + blockIdx.z*blockDim.z - input.nguard;
#else
    int k = 0;
#endif
    double xyz[3];
    double dx[3];
    dx[0] = (input.bounds[1] - input.bounds[0])/input.nxb[0];
    dx[1] = (input.bounds[3] - input.bounds[2])/input.nxb[1];
#if(IS3D)
    dx[2] = (input.bounds[5] - input.bounds[4])/input.nxb[2];
#else
    dx[2] = 0.0;
#endif

    if (i < (input.nxb[0]+input.nguard) && j < (input.nxb[1]+input.nguard) && k < (input.nxb[2]+input.nguard))
    {
        xyz[0] = input.bounds[0] + (i - input.nguard + 0.5)*dx[0];
        xyz[1] = input.bounds[2] + (j - input.nguard + 0.5)*dx[1];
#if(IS3D)
        xyz[2] = input.bounds[4] + (k - input.nguard + 0.5)*dx[2];
#else
        xyz[2] = 0.0;
#endif
        flow[bidx(0, i, j, k, lb, input)] = 0.0;
        flow[bidx(1, i, j, k, lb, input)] = 0.0;
        flow[bidx(2, i, j, k, lb, input)] = 0.0;
        flow[bidx(3, i, j, k, lb, input)] = 0.0;
#if(IS3D)
        flow[bidx(4, i, j, k, lb, input)] = 0.0;
#endif
    }
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
    gridConf.x = (numcells[0] + BLOCK_SIZE - 1)/BLOCK_SIZE;
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
    
    for (int lb = 0; lb < input.lnblocks; lb++)
    {
        K_TestFunctionsGpu<<<gridConf, blockConf>>>(flow, input, lb);
        CuCheck(cudaPeekAtLastError());
    }
    CuCheck(cudaDeviceSynchronize());
}