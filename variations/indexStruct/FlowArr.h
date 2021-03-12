#ifndef HIDX_FLOW_ARR_H
#define HIDX_FLOW_ARR_H
#include "InputClass.h"
#include <cuda.h>
#include "cuda_device_runtime_api.h"
#include "cuda_runtime_api.h"
#include <cuda_runtime.h>
#include "CuErr.h"
struct FlowArr
{
    double* data;
    int nv, ni, nj, nk;
    int di, dj, dk;
    bool isCpu;
    FlowArr(int nvIn, int niIn, int njIn, int nkIn, int nlbIn, int nguardI, int nguardJ, int nguardK, device::device dev)
    {
        di = nguardI;
        dj = nguardJ;
        dk = nguardK;
        isCpu = (dev==device::cpu);
        size_t totalSize = sizeof(double)*(nvIn)*(niIn+2*nguardI)*(njIn+2*nguardJ)*(nkIn+2*nguardK)*(nlbIn);
        if (isCpu)
        {
            data = (double*)malloc(totalSize);
        }
        else
        {
            CuCheck(cudaMalloc((void**)(&data), totalSize));
        }
        nv = nvIn;
        ni = niIn + 2*nguardI;
        nj = njIn + 2*nguardJ;
        nk = nkIn + 2*nguardK;
    }
    
    double& operator () (int v, int i, int j, int k, int lb)
    {
        return *(data + v + (i+di)*nv + (j+dj)*nv*ni + (k+dk)*nv*ni*nj + lb*nv*ni*nj*nk);
    }
    
    ~FlowArr(void)
    {
        if (isCpu)
        {
            free(data);
        }
        else
        {
            CuCheck(cudaFree(data));
        }
    }
};

#endif