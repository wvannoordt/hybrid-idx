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
#include "Result.h"
#include <chrono>
#include "FlowArr.h"

bool HasEnding(std::string const &fullString, std::string const &ending)
{
    if (fullString.length() >= ending.length())
    {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    }
    else
    {
        return false;
    }
}

#ifndef OPTL
#define OPTL 0
#endif

InputClass input;

int main(int argc, char** argv)
{
    std::string inputfile = "input.ptl";
    
    for (int i = 0; i < argc; i++)
    {
        std::string argStr(argv[i]);
        if (HasEnding(argStr, ".ptl"))
        {
            inputfile = argStr;
        }
    }
    
    // MPI_Init(&argc, &argv);
    int mypeno, noprocs;
    // MPI_Comm_rank(MPI_COMM_WORLD, &mypeno);
    // MPI_Comm_size(MPI_COMM_WORLD, &noprocs);
    int* nxbTmp;
    double* boundsTmp;
    std::string outfile;
    PTL::PropertyTree ptree;
    PTL::Interactive i(argc, argv, &ptree);
    ptree["centOrder"].MapTo(&input.centOrder) = new PTL::PTLInteger(2, "order of central scheme");
    ptree["nxb"].MapTo(&nxbTmp) = new PTL::PTLStaticIntegerArray(DIM, "dimensions of blocks", [](int i){return 16;});
    ptree["lnblocks"].MapTo(&input.lnblocks) = new PTL::PTLInteger(10, "total number of blocks");
    ptree["numSteps"].MapTo(&input.numSteps) = new PTL::PTLInteger(1000, "number of timesteps");
    ptree["nguard"].MapTo(&input.nguard) = new PTL::PTLInteger(4, "number of guard cells");
    ptree["bounds"].MapTo(&boundsTmp) = new PTL::PTLStaticDoubleArray(2*DIM, "bounds of the block", [](int i){return (double)(0 + (i%2));});
    ptree["Rgas"].MapTo(&input.Rgas)  = new PTL::PTLDouble(287.0, "Gas constant");
    ptree["gamma"].MapTo(&input.gamma)  = new PTL::PTLDouble(1.4, "Specific heat ratio");
    ptree["device"].MapTo(&input.dev)  = new PTL::PTLAutoEnum(device::cpu, deviceStr, "The device to run on");
    ptree["outputGuards"].MapTo(&input.outputGuards) = new PTL::PTLBoolean(false, "Output the guards");
    ptree["outputError"].MapTo(&input.outputError) = new PTL::PTLBoolean(false, "Output the error file");
    ptree["resultFile"].MapTo(&outfile) = new PTL::PTLString("output/default.json", "Output the error file");
    mypeno = 0;
    bool isGpu = (input.dev == device::gpu);
    ptree.Read(inputfile);
    ptree.StrictParse();
    mypenoG = mypeno;
    hasPrintedGp = false;
    for (int i = 0; i < DIM; i++)
    {
        input.bounds[2*i] = boundsTmp[2*i];
        input.bounds[2*i+1] = boundsTmp[2*i+1];
        input.nxb[i] = nxbTmp[i];
    }
    // std::cout << "Hello from " << mypeno << " of " << noprocs << std::endl;
    // MPI_Barrier(MPI_COMM_WORLD);
    
    size_t blockSize = 1;
    for (int i = 0; i < DIM; i++) blockSize*=(2*input.nguard + input.nxb[i]);
    size_t totalSize = sizeof(double) * blockSize * (2 + DIM) * input.lnblocks;
    if (mypeno == 0) std::cout << totalSize << " (bytes)" << std::endl;
    if (mypeno == 0) std::cout << totalSize/sizeof(double) << " (elements)" << std::endl;
    if (mypeno == 0) std::cout << totalSize/(sizeof(double)*(2+DIM)) << " (cells)" << std::endl;
    FlowArr cpuFlow(2+DIM, input.nxb[0], input.nxb[1], 1-IS3D+IS3D*input.nxb[DIM-1], input.lnblocks, input.nguard, input.nguard, input.nguard, device::cpu);
    cpuFlow.Alloc(input.lnblocks);
    FlowArr gpuFlow(2+DIM, input.nxb[0], input.nxb[1], 1-IS3D+IS3D*input.nxb[DIM-1], input.lnblocks, input.nguard, input.nguard, input.nguard, device::gpu);
    if (isGpu)
    {
        gpuFlow.Alloc(input.lnblocks);
        CuCheck(cudaMalloc((void**)(&gpuFlow.data), gpuFlow.totalSize));
    }
    
    FlowArr cpuErr(2+DIM, input.nxb[0], input.nxb[1], 1-IS3D+IS3D*input.nxb[DIM-1], input.lnblocks, input.nguard, input.nguard, input.nguard, device::cpu);
    cpuErr.Alloc(input.lnblocks);
    FlowArr gpuErr(2+DIM, input.nxb[0], input.nxb[1], 1-IS3D+IS3D*input.nxb[DIM-1], input.lnblocks, input.nguard, input.nguard, input.nguard, device::gpu);
    FlowArr gpuErrMirror(2+DIM, input.nxb[0], input.nxb[1], 1-IS3D+IS3D*input.nxb[DIM-1], input.lnblocks, input.nguard, input.nguard, input.nguard, device::cpu);
    if (isGpu)
    {
        gpuErr.Alloc(input.lnblocks);
        CuCheck(cudaMalloc((void**)(&gpuErr.data), gpuErr.totalSize));
        gpuErrMirror.Alloc(input.lnblocks);
    }
    
    InitCpu(cpuFlow, cpuErr, input);
    if (isGpu) InitGpu(gpuFlow, gpuErr, input);
    int numZer = std::to_string(input.numSteps-1).length();
    double elapsedTime = 0.0;
    for (int nt = 0; nt < input.numSteps; nt++)
    {
        auto start = std::chrono::high_resolution_clock::now();
        if (input.dev == device::gpu) {ConvGpu(gpuFlow, gpuErr, input);}
        else {ConvCpu(cpuFlow, cpuErr, input);}
        // MPI_Barrier(MPI_COMM_WORLD);
        
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        double timeMS = 1000*elapsed.count();
        elapsedTime += timeMS;
        if (mypenoG == 0) std::cout << ("nt: " + zfill(nt, numZer) + "/" + zfill(input.numSteps-1, numZer) + " time: " + std::to_string(timeMS) + " ms") << std::endl;
    }
    bool pass = false;
    // MPI_Barrier(MPI_COMM_WORLD);
    if (mypenoG == 0)
    {
        std::cout << "Average timestep: " + std::to_string(elapsedTime/input.numSteps) + " ms" << std::endl;
    }
    if ((input.dev == device::cpu) && (mypenoG==0) && input.outputError)
    {
        pass = Output(cpuErr, input, 0, "output/cpu.vtk");
    }
    if ((input.dev == device::gpu) && (mypenoG==0) && input.outputError)
    {
        GCopy(gpuErrMirror, gpuErr, totalSize);
        pass = Output(gpuErrMirror, input, 0, "output/gpu.vtk");
    }
    if (mypeno == 0) std::cout << "Cleaning up" << std::endl;
    
    if (mypenoG==0)
    {
        Result result;
        result.avgStepTime = elapsedTime/input.numSteps;
        std::ofstream myfile;
        myfile.open(outfile.c_str());
        myfile << "{\n";
        myfile << "    \"optLevel\": " << ((OPTL>0)?("\"opt\""):("\"dbg\"")) << ",\n";
        myfile << "    \"cpuKernel\": \"" << GetCpuKernelDescription() << "\",\n";
        myfile << "    \"gpuKernel\": \"" << GetGpuKernelDescription() << "\",\n";
        myfile << "    \"pass\": " << ((pass)?("true"):("false")) << ",\n";
        input.WriteJson(myfile);
        result.WriteJson(myfile);
        myfile << "}";
        myfile.close();
    }
    
    // MPI_Finalize();
    return pass?0:187;
}
