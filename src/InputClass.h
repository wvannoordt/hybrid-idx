#ifndef INPUT_CLS_H
#define INPUT_CLS_H
#include "Config.h"
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
};
#endif
