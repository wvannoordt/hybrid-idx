#ifndef GPU_KERNEL_H
#define GPU_KERNEL_H
#include "InputClass.h"
#include "Config.h"
#include "Glob.h"
struct Coef_t
{
    double c[4];
};
void InitGpu(double* flow, double* err, const InputClass& input);
void ConvGpu(double* flow, double* err, const InputClass& input);
#endif