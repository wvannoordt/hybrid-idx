#ifndef GPU_KERNEL_H
#define GPU_KERNEL_H
#include "InputClass.h"
#include "Config.h"
#include "Glob.h"
#include <string>
#include "FlowArr.h"
struct Coef_t
{
    double c[4];
};
void InitGpu(FlowArr& flow, FlowArr& err, const InputClass& input);
void ConvGpu(FlowArr& flow, FlowArr& err, const InputClass& input);
void GCopy(FlowArr& cTarget, FlowArr& gTarget, size_t size);
std::string GetGpuKernelDescription(void);
#endif