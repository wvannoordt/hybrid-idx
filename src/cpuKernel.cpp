#include "cpuKernel.h"
#include "Idx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
void InitCpu(double* flow, const InputClass& input)
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
                    for (int d = 0; d < DIM; d++) xyz[d] = input.bounds[2*d] + (ijk[d] + 0.5) * dx[d];
                    flow[bidx(0, i, j, k, lb, input)] = xyz[0];
                    flow[bidx(1, i, j, k, lb, input)] = xyz[1];
                    flow[bidx(2, i, j, k, lb, input)] = xyz[2];
                    flow[bidx(3, i, j, k, lb, input)] = xyz[0];
#if(IS3D)
                    flow[bidx(4, i, j, k, lb, input)] = xyz[1];
#endif
                }
            }
        }
    }
}

void ConvCpu(double* flow, const InputClass& input)
{
    
    int imin = 0;
    int imax = input.nxb[0];
    int jmin = 0;
    int jmax = input.nxb[1];
    int kmin = 0;
    int kmax = IS3D*input.nxb[1+IS3D] + (1-IS3D);
    for (int k = kmin; k < kmax; k++)
    {
        for (int j = jmin; j < jmax; j++)
        {
            for (int i = imin; i < imax; i++)
            {
                
            }
        }
    }
}

void OutputCpu(double* flow, const InputClass& input, int lb)
{
    int imin = -input.nguard;
    int imax = (input.nxb[0] + input.nguard)+1;
    int jmin = -input.nguard;
    int jmax = (input.nxb[1] + input.nguard)+1;
    int kmin = IS3D*(-input.nguard);
    int kmax = IS3D*(input.nxb[1+IS3D] + input.nguard) + (1-IS3D) + IS3D;
    double xyz[3];
    xyz[0] = 0.0; xyz[1] = 0.0; xyz[2] = 0.0;
    double dx[DIM];
    int ijk[3];
    for (int d = 0; d < DIM; d++) dx[d] = (input.bounds[2*d+1] - input.bounds[2*d])/input.nxb[d];
    
    std::string filename = "output/block" + zfill(lb, 4) + "p" + zfill(mypenoG, 3) + ".vtk";
    std::ofstream myfile;
    myfile.open(filename.c_str());
    myfile << "# vtk DataFile Version 3.0" << std::endl;
    myfile << "GME to the moon" << std::endl;
    myfile << "ASCII" << std::endl;
    myfile << "DATASET STRUCTURED_GRID" << std::endl;
    myfile << "DIMENSIONS " << imax-imin << " " << jmax-jmin << " " << kmax-kmin << std::endl;
    myfile << "POINTS " << (imax-imin)*(jmax-jmin)*(kmax-kmin) << " double" << std::endl;
    for (int k = kmin; k < kmax; k++)
    {
        ijk[2] = k;
        for (int j = jmin; j < jmax; j++)
        {
            ijk[1] = j;
            for (int i = imin; i < imax; i++)
            {
                ijk[0] = i;
                for (int d = 0; d < DIM; d++) xyz[d] = input.bounds[2*d] + (ijk[d]) * dx[d];
                myfile << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
            }
        }
    }
    imax = (input.nxb[0] + input.nguard);
    jmax = (input.nxb[1] + input.nguard);
    kmax = IS3D*(input.nxb[1+IS3D] + input.nguard) + (1-IS3D);
    myfile << "CELL_DATA " << (imax-imin)*(jmax-jmin)*(kmax-kmin) << std::endl;
    int v = 0;
    std::vector<std::string> names;
    names.push_back("P");
    names.push_back("T");
    names.push_back("U");
    names.push_back("V");
#if(IS3D)
    names.push_back("W");
#endif
    for (const auto s: names)
    {
        myfile << "SCALARS " << s << " double" << std::endl;
        myfile << "LOOKUP_TABLE default" << std::endl;
        for (int k = kmin; k < kmax; k++)
        {
            ijk[2] = k;
            for (int j = jmin; j < jmax; j++)
            {
                ijk[1] = j;
                for (int i = imin; i < imax; i++)
                {
                    ijk[0] = i;
                    myfile << flow[bidx(v, i, j, k, lb, input)] << std::endl;
                }
            }
        }
        v++;
    }
    myfile.close();
}