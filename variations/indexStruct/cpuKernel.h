#ifndef CPU_KERNEL_H
#define CPU_KERNEL_H
#include "InputClass.h"
#include "Config.h"
#include "Glob.h"
#include <string>
#include "FlowArr.h"
void InitCpu(FlowArr& flow, FlowArr& err, const InputClass& input);
bool Output(FlowArr& flow, const InputClass& input, int lb, std::string filename);
void ConvCpu(FlowArr& flow, FlowArr& err, const InputClass& input);
std::string GetCpuKernelDescription(void);
#endif