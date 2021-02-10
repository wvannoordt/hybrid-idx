#ifndef CPU_KERNEL_H
#define CPU_KERNEL_H
#include "InputClass.h"
#include "Config.h"
#include "Glob.h"
#include <string>
void InitCpu(double* flow, double* err, const InputClass& input);
void Output(double* flow, const InputClass& input, int lb, std::string filename);
void ConvCpu(double* flow, double* err, const InputClass& input);
#endif