#ifndef GPU_KERNEL_H
#define GPU_KERNEL_H
#include "InputClass.h"
#include "Config.h"
#include "Glob.h"
#include <string>
struct Coef_t
{
    double c[4];
};
void InitGpu(double* flow, double* err, const InputClass& input);
void ConvGpu(double* flow, double* err, const InputClass& input);
void GCopy(double* cTarget, double* gTarget, size_t size);
std::string GetGpuKernelDescription(void);
#endif