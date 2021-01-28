#include "cpuKernel.h"
#include <iostream>
void TestFunctionsCpu(double* flow, const InputClass& input)
{
    int imin = 0;
    int imax = (input.nxb[0] + 2*input.nguard);
    int jmin = 0;
    int jmax = (input.nxb[1] + 2*input.nguard);
    int kmin = 0;
    int kmax = IS3D*(input.nxb[1+IS3D] + 2*input.nguard);
    int num = 0;
    
    double xyz[3];
    xyz[0] = 0.0; xyz[1] = 0.0; xyz[2] = 0.0;
    double dx[DIM];
    int ijk[3];
    for (int d = 0; d < DIM; d++) dx[d] = (input.bounds[2*d+1] - input.bounds[2*d])/input.nxb[d];
    for (int k = kmin; k < kmax; k++)
    {
        ijk[2] = k;
        for (int j = jmin; j < jmax; j++)
        {
            ijk[1] = j;
            for (int i = imin; i < imax; i++)
            {
                ijk[0] = i;
                for (int d = 0; d < DIM; d++) xyz[d] = input.bounds[2*d] + (ijk[d] - input.nguard + 0.5) * dx[d];
                
                // do stuff with xyz
            }
        }
    }
}