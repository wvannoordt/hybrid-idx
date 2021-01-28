#include "cpuKernel.h"
#include "Idx.h"
#include <iostream>
void TestFunctionsCpu(double* flow, const InputClass& input)
{
    int imin = -input.nguard;
    int imax = (input.nxb[0] + input.nguard);
    int jmin = -input.nguard;
    int jmax = (input.nxb[1] + input.nguard);
    int kmin = IS3D*(-input.nguard);
    int kmax = IS3D*(input.nxb[1+IS3D] + input.nguard) + (1-IS3D);
    int num = 0;
    
    double xyz[3];
    xyz[0] = 0.0; xyz[1] = 0.0; xyz[2] = 0.0;
    double dx[DIM];
    int ijk[3];
    for (int d = 0; d < DIM; d++) dx[d] = (input.bounds[2*d+1] - input.bounds[2*d])/input.nxb[d];
    for (int lb = 0; lb < input.lnblocks; lb++)
    {
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
                    flow[bidx(0, i, j, k, lb, input)] = 0.0;
                    flow[bidx(1, i, j, k, lb, input)] = 0.0;
                    flow[bidx(2, i, j, k, lb, input)] = 0.0;
                    flow[bidx(3, i, j, k, lb, input)] = 0.0;
                    flow[bidx(4, i, j, k, lb, input)] = 0.0;
                }
            }
        }
    }
}