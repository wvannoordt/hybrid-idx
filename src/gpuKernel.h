#ifndef GPU_KERNEL_H
#define GPU_KERNEL_H
#include "InputClass.h"
#include "Config.h"
#include "Glob.h"
void InitGpu(double* flow, const InputClass& input);
void ConvGpu(double* flow, const InputClass& input, int lb);
#endif