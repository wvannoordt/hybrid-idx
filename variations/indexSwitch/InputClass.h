#ifndef INPUT_CLS_H
#define INPUT_CLS_H
#include "Config.h"
#include <string>
#include <fstream>
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
    
    void WriteJson(std::ofstream& fs)
    {
        fs << "    \"centOrder\": " << centOrder << ",\n";
        fs << "    \"lnblocks\": " << lnblocks << ",\n";
        fs << "    \"nguard\": " << nguard << ",\n";
        fs << "    \"numSteps\": " << numSteps << ",\n";
        fs << "    \"Rgas\": " << Rgas << ",\n";
        fs << "    \"gamma\": " << gamma << ",\n";
        fs << "    \"dev\": " << dev << ",\n";
#if(IS3D)
        fs << "    \"nxb\": [" << nxb[0] << ", " << nxb[1] << ", " << nxb[2] << "],\n";
#else
        fs << "    \"nxb\": [" << nxb[0] << ", " << nxb[1] << "],\n";
#endif
    }
};
#endif
