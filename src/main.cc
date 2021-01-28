#include "PTL.h"
#include "mpi.h"
#include "Config.h"
#include "InputClass.h"
#include "cpuKernel.h"
#include "gpuKernel.h"
#include <cuda.h>
#include "cuda_device_runtime_api.h"
#include "cuda_runtime_api.h"
#include <cuda_runtime.h>

InputClass input;

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mypeno, noprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &mypeno);
    MPI_Comm_size(MPI_COMM_WORLD, &noprocs);
    int* nxbTmp;
    double* boundsTmp;
    PTL::PropertyTree ptree;
    ptree["centOrder"].MapTo(&input.centOrder) = new PTL::PTLInteger(2, "order of central scheme");
    ptree["nxb"].MapTo(&nxbTmp) = new PTL::PTLStaticIntegerArray(DIM, "dimensions of blocks", [](int i){return 16;});
    ptree["lnblocks"].MapTo(&input.lnblocks) = new PTL::PTLInteger(10, "total number of blocks");
    ptree["numSteps"].MapTo(&input.numSteps) = new PTL::PTLInteger(1000, "number of timesteps");
    ptree["nguard"].MapTo(&input.nguard) = new PTL::PTLInteger(4, "number of guard cells");
    ptree["bounds"].MapTo(&boundsTmp) = new PTL::PTLStaticDoubleArray(2*DIM, "bounds of the block", [](int i){return (double)(0 + (i%2));});
    
    ptree.Read("input.ptl");
    ptree.StrictParse();
    for (int i = 0; i < DIM; i++)
    {
        input.bounds[2*i] = boundsTmp[2*i];
        input.bounds[2*i+1] = boundsTmp[2*i+1];
        input.nxb[i] = nxbTmp[i];
    }
    std::cout << "Hello from " << mypeno << " of " << noprocs << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    
    size_t blockSize = 1;
    for (int i = 0; i < DIM; i++) blockSize*=(2*input.nguard + input.nxb[i]);
    size_t totalSize = sizeof(double) * blockSize * (2 + DIM);
    if (mypeno == 0) std::cout << totalSize << " (bytes)" << std::endl;
    double* cpuFlow = (double*)malloc(totalSize);
    double* gpuFlow;
    cudaMalloc(&gpuFlow, totalSize);
    
    TestFunctionsCpu(cpuFlow, input);
    TestFunctionsCpu(gpuFlow, input);
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (mypeno == 0) std::cout << "Cleaning up" << std::endl;
    free (cpuFlow);
    cudaFree(gpuFlow);
    MPI_Finalize();
    
}