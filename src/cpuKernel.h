#ifndef CPU_KERNEL_H
#define CPU_KERNEL_H
#include "InputClass.h"
#include "Config.h"
#include "Glob.h"
void InitCpu(double* flow, double* err, const InputClass& input);
void OutputCpu(double* flow, const InputClass& input, int lb);
void ConvCpu(double* flow, double* err, const InputClass& input);
#endif