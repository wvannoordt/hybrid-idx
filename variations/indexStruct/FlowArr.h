#ifndef HIDX_FLOW_ARR_H
#define HIDX_FLOW_ARR_H
#include "InputClass.h"
#include <cuda.h>
#include "cuda_device_runtime_api.h"
#include "cuda_runtime_api.h"
#include <cuda_runtime.h>
#include "CuErr.h"

#ifdef __CUDACC__
#define __share __device__ __host__
#else
#define __share
#endif

struct FlowArr
{
    double* data;
    int nv, ni, nj, nk;
    int di, dj, dk;
    bool isCpu;
    size_t totalSize;
    FlowArr(int nvIn, int niIn, int njIn, int nkIn, int nlbIn, int nguardI, int nguardJ, int nguardK, device::device dev)
    {
        di = nguardI;
        dj = nguardJ;
        dk = nguardK;
        isCpu = (dev==device::cpu);
        nv = nvIn;
        ni = niIn + 2*nguardI;
        nj = njIn + 2*nguardJ;
        nk = nkIn + 2*nguardK;
        totalSize = sizeof(double)*(nv)*(ni)*(nj)*(nk)*(nlbIn);
    }
    
    void Alloc(int lbin)
    {
        if (isCpu)
        {
            data = (double*)malloc(totalSize);
        }
    }
    
    __share double& operator () (int v, int i, int j, int k, int lb)
    {
        return *(data + v + (i+di)*nv + (j+dj)*nv*ni + (k+dk)*nv*ni*nj + lb*nv*ni*nj*nk);
    }
    
    ~FlowArr(void)
    {
        std::cout << "close" << std::endl;
    }
};

#endif