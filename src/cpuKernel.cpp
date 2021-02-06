#include "cpuKernel.h"
#include "Idx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
void InitCpu(double* flow, double* err, const InputClass& input)
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
                    err[bidx(0, i, j, k, lb, input)] = 0.0;
                }
            }
        }
    }
}

#define stencilIdx(v,j) ((v)+(5+DIM)*(j))

#define f_DivSplit(q,j,l,v1)         (0.500*(q[stencilIdx((v1),(j))] + q[stencilIdx((v1),(j)+(l))]))
#define fg_QuadSplit(q,j,l,v1,v2)    (0.250*(q[stencilIdx((v1),(j))] + q[stencilIdx((v1),(j)+(l))])*(q[stencilIdx((v2),(j))] + q[stencilIdx((v2),(j)+(l))]))
#define fg_CubeSplit(q,j,l,v1,v2,v3) (0.125*(q[stencilIdx((v1),(j))] + q[stencilIdx((v1),(j)+(l))])*(q[stencilIdx((v2),(j))] + q[stencilIdx((v2),(j)+(l))])*(q[stencilIdx((v3),(j))] + q[stencilIdx((v3),(j)+(l))]))
#define fg_DivSplit(q,j,l,v1,v2)     (0.500*(q[stencilIdx((v1),(j)+(l))]*q[stencilIdx((v2),(j))]) + (q[stencilIdx((v1),(j))]*q[stencilIdx((v2),(j)+(l))]))

void ConvCpu(double* flow, double* err, const InputClass& input)
{
    double centerCoef[4] = {0.0};
    switch (input.centOrder)
    {
        case 2: {centerCoef[0] = 1.0/2.0; break;}
        case 4: {centerCoef[0] = 2.0/3.0; centerCoef[1] = -1.0/12.0; break;}
        case 6: {centerCoef[0] = 3.0/4.0; centerCoef[1] = -3.0/20.0; centerCoef[2] = 1.0/60.0 ;break;}
        case 8: {centerCoef[0] = 4.0/5.0; centerCoef[1] = -1.0/5.0 ; centerCoef[2] = 4.0/105; centerCoef[3] = -1.0/280.0; break;}
        default: {std::cout << "Bad central scheme order." << std::endl; abort();}
    }

    // p, T, u, v, w
    int imin = 0;
    int imax = input.nxb[0];
    int jmin = 0;
    int jmax = input.nxb[1];
    int kmin = 0;
    int kmax = IS3D*input.nxb[1+IS3D] + (1-IS3D);
    int stencilWid = input.centOrder/2;
    int dijk[3];

    double invdx[DIM] = {0};

    for (int lb = 0; lb < input.lnblocks; lb++)
    {        
        for (int i = 0; i < DIM; i++)
        {
            invdx[i] = (double) (input.nxb[i])/(input.bounds[2*i+1] - input.bounds[2*i]);
        }

        for (int idir = 0; idir < DIM; idir ++)
        {
            for (int d = 0; d < DIM; d++) dijk[d] = 0;
            dijk[idir] = 1;
            for (int k = kmin; k < kmax; k++)
            {
                for (int j = jmin; j < jmax; j++)
                {
                    for (int i = imin; i < imax; i++)
                    {
                        double stencilData[9*(5+DIM)]; //ie,ke,rho,P,T,u,v,w
                        double rhs[2+DIM] = {0.0};
                        // fluxes
                        double C[2]     = {0.0};
                        double M[DIM*2] = {0.0};
                        double PGRAD[2] = {0.0};
                        double KE[2]    = {0.0};
                        double IE[2]    = {0.0};
                        double PDIFF[2] = {0.0};

                        for (int n = 0; n < input.centOrder + 1; n++)
                        {
                            for (int v = 3; v < (5+DIM); v++)
                            {
                                int ii = i+dijk[0]*(n-stencilWid);
                                int jj = j+dijk[1]*(n-stencilWid);
                                int kk = k+dijk[2]*(n-stencilWid);
                                int h = bidx(v-3, ii, jj, kk, lb, input);
                                stencilData[stencilIdx(v,n)] = flow[h];
                            }
                            // IE = P/(rho*(gamma - 1))
                            stencilData[stencilIdx(0,n)] = stencilData[stencilIdx(3,n)]/(stencilData[stencilIdx(2,n)]*(input.gamma - 1.0));
                            stencilData[stencilIdx(1,n)] = 0.0;

                            // Not needed per se starts
                            for (int vel_comp = 0; vel_comp < DIM; vel_comp ++)
                            {
                                stencilData[stencilIdx(1,n)] += 0.5*stencilData[stencilIdx(5+vel_comp,n)]*stencilData[stencilIdx(5+vel_comp,n)];
                            }
                            // Not needed per se ends

                            stencilData[stencilIdx(2,n)] = stencilData[stencilIdx(3,n)]/(input.Rgas*stencilData[stencilIdx(4,n)]);
                        }
                        // Mass conservation                              
                        for (int l = 1; l <= stencilWid; l++)
                        {
                            double al = centerCoef[l-1];
                            int jf = stencilWid;
                            for (int m = 0; m <= (l-1); m++)
                            {
                                C[1] += 2.0*centerCoef[l-1]*fg_QuadSplit(stencilData,jf-m, l,2,5+idir);
                                C[0] += 2.0*centerCoef[l-1]*fg_QuadSplit(stencilData,jf+m,-l,2,5+idir);
                                for (int idir_mom = 0; idir_mom < DIM; idir_mom++)
                                {
                                    M[idir_mom      ] = 2.0*centerCoef[l-1]*fg_CubeSplit(stencilData,jf-m, l,3,5+idir,5+idir_mom);
                                    M[idir_mom + DIM] = 2.0*centerCoef[l-1]*fg_CubeSplit(stencilData,jf+m,-l,3,5+idir,5+idir_mom);
                                }

                                PGRAD[1] += 2.0*centerCoef[l-1]*f_DivSplit(stencilData,jf-m, l,3);
                                PGRAD[0] += 2.0*centerCoef[l-1]*f_DivSplit(stencilData,jf+m,-l,3);

                                for (int vel_comp = 0;  vel_comp < DIM; vel_comp ++)
                                {
                                    KE[1] += 2.0*centerCoef[l-1]*fg_QuadSplit(stencilData,jf-m, l,2,5+idir)*0.5*(stencilData[stencilIdx(5+vel_comp,jf-m)]*stencilData[stencilIdx(5+vel_comp,jf-m+l)]);
                                    KE[0] += 2.0*centerCoef[l-1]*fg_QuadSplit(stencilData,jf+m,-l,2,5+idir)*0.5*(stencilData[stencilIdx(5+vel_comp,jf+m)]*stencilData[stencilIdx(5+vel_comp,jf+m-l)]);
                                }

                                IE[1] += 2.0*centerCoef[l-1]*fg_CubeSplit(stencilData,jf-m, l,3,0,5+idir);
                                IE[0] += 2.0*centerCoef[l-1]*fg_CubeSplit(stencilData,jf+m,-l,3,0,5+idir);

                                PDIFF[1] += 2.0*centerCoef[l-1]*fg_DivSplit(stencilData,jf-m, l,5+idir,3);
                                PDIFF[0] += 2.0*centerCoef[l-1]*fg_DivSplit(stencilData,jf-m, l,5+idir,3);
                            }
                        }

                        rhs[0] += -invdx[idir]*(C[1] - C[0]);
                        rhs[1] += -invdx[idir]*(IE[1] + KE[1] + PDIFF[1] - IE[0] - KE[0] - PDIFF[0]);
                        rhs[2+idir] += -invdx[idir]*(PGRAD[1] - PGRAD[0]);
                        for (int rhs_vel_comp = 0; rhs_vel_comp < DIM; rhs_vel_comp++)
                        {
                            rhs[2+rhs_vel_comp] += -invdx[idir]*(M[1] - M[0]);
                        }
                        for (int y = 0; y < 2+DIM; y++) err[bidx(0, i, j, k, lb, input)] += rhs[y];
                    }
                }
            }
            dijk[idir] = 0;
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
