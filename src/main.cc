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
#include "Glob.h"
#include "CuErr.h"
#include <chrono>



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
    ptree["Rgas"].MapTo(&input.Rgas)  = new PTL::PTLDouble(287.0, "Gas constant");
    ptree["gamma"].MapTo(&input.Rgas)  = new PTL::PTLDouble(1.4, "Gamma (Specific heat ratio)");
    ptree["device"].MapTo(&input.dev)  = new PTL::PTLAutoEnum(device::cpu, deviceStr, "The device to run on");
    ptree["outputGuards"].MapTo(&input.outputGuards) = new PTL::PTLBoolean(false, "Output the guards");
    
    ptree.Read("input.ptl");
    ptree.StrictParse();
    mypenoG = mypeno;
    hasPrintedGp = false;
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
    size_t totalSize = sizeof(double) * blockSize * (2 + DIM) * input.lnblocks;
    if (mypeno == 0) std::cout << totalSize << " (bytes)" << std::endl;
    if (mypeno == 0) std::cout << totalSize/sizeof(double) << " (elements)" << std::endl;
    double* cpuFlow = (double*)malloc(totalSize);
    double* gpuFlow = 0;
    CuCheck(cudaMalloc((void**)(&gpuFlow), totalSize));
    
    double* cpuErr = (double*)malloc(totalSize);
    double* gpuErr = 0;
    CuCheck(cudaMalloc((void**)(&gpuErr), totalSize));
    
    InitCpu(cpuFlow, cpuErr, input);
    InitGpu(gpuFlow, gpuErr, input);
    int numZer = std::to_string(input.numSteps-1).length();
    double elapsedTime = 0.0;
    for (int nt = 0; nt < input.numSteps; nt++)
    {
        auto start = std::chrono::high_resolution_clock::now();
        if (input.dev == device::gpu) {ConvGpu(gpuFlow, gpuErr, input);}
        else {ConvCpu(cpuFlow, cpuErr, input);}
        MPI_Barrier(MPI_COMM_WORLD);
        
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        double timeMS = 1000*elapsed.count();
        elapsedTime += timeMS;
        if (mypenoG == 0) std::cout << ("nt: " + zfill(nt, numZer) + "/" + zfill(input.numSteps-1, numZer) + " time: " + std::to_string(timeMS) + " ms") << std::endl;
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (mypenoG == 0)
    {
        std::cout << "Average timestep: " + std::to_string(elapsedTime/input.numSteps) + " ms" << std::endl;
    }
    if ((input.dev == device::cpu) && (mypenoG==0))
    {
        OutputCpu(cpuErr, input, 0);
    }
    if ((input.dev == device::gpu) && (mypenoG==0))
    {
        // OutputGpu(gpuErr, input, 0);
    }
    if (mypeno == 0) std::cout << "Cleaning up" << std::endl;
    free(cpuFlow);
    CuCheck(cudaFree(gpuFlow));
    free(cpuErr);
    CuCheck(cudaFree(gpuErr));
    MPI_Finalize();
    
}
