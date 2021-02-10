#ifndef INPUT_CLS_H
#define INPUT_CLS_H
#include "Config.h"
#include <string>
#include "PTL.h"
namespace device{enum device{cpu,gpu};}
inline static std::string deviceStr(int d)
{
    switch (d)
    {
        case device::cpu: return "cpu";
        case device::gpu: return "gpu";
    }
    return PTL_AUTO_ENUM_TERMINATOR;
}

struct InputClass
{
    int centOrder;
    int nxb[DIM];
    int lnblocks;
    int nguard;
    int numSteps;
    double bounds[2*DIM];
    double Rgas;
    double gamma;
    int dev;
    bool outputGuards;
    bool outputError;
};
#endif
