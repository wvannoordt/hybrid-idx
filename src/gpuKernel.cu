#include "gpuKernel.h"
#include "Idx.h"
#include "CuErr.h"
#include <iostream>
__global__ void K_Init(double* flow, double* err, const InputClass input, const int lb)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x - input.nguard;
    int j = threadIdx.y + blockIdx.y*blockDim.y - input.nguard;
#if(IS3D)    
    int k = threadIdx.z + blockIdx.z*blockDim.z - input.nguard;
#else
    int k = 0;
#endif
    double xyz[3];
    double dx[3];
    dx[0] = (input.bounds[1] - input.bounds[0])/input.nxb[0];
    dx[1] = (input.bounds[3] - input.bounds[2])/input.nxb[1];
#if(IS3D)
    dx[2] = (input.bounds[5] - input.bounds[4])/input.nxb[2];
#else
    dx[2] = 0.0;
#endif

    if (i < (input.nxb[0]+input.nguard) && j < (input.nxb[1]+input.nguard) && k < (input.nxb[2]+input.nguard))
    {
        xyz[0] = input.bounds[0] + (i + 0.5)*dx[0];
        xyz[1] = input.bounds[2] + (j + 0.5)*dx[1];
#if(IS3D)
        xyz[2] = input.bounds[4] + (k + 0.5)*dx[2];
#else
        xyz[2] = 0.0;
#endif
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

void InitGpu(double* flow, double* err, const InputClass& input)
{
    //Global memory:            7.907288 GB
    //Shared memory per block:  48.000000 KB
    //Warp size:                32￼
    //Max threads per block:    1024
    //Max thread dimension:     1,024  x  1,024  x  64
    //Max grid size:            2,147,483,647  x  65,535  x  65,535
    //Total constant memory:    64.000000 KB
    
    dim3 blockConf;
    blockConf.x = BLOCK_SIZE;
    blockConf.y = BLOCK_SIZE;
#if(IS3D)
    blockConf.z = BLOCK_SIZE;
#endif
    dim3 gridConf;
    int numcells[DIM];
    for (int i = 0; i < DIM; i++) {numcells[i] = input.nxb[i] + 2*input.nguard;}
    gridConf.x = (numcells[0] + BLOCK_SIZE - 1)/BLOCK_SIZE;
    gridConf.y = (numcells[1] + BLOCK_SIZE - 1)/BLOCK_SIZE;
#if(IS3D)
    gridConf.z = (numcells[2] + BLOCK_SIZE - 1)/BLOCK_SIZE;
#endif

    if (mypenoG==0 && !hasPrintedGp)
    {
        hasPrintedGp = true;
        std::cout << "GP Config:\nblock: " << blockConf.x << " x " << blockConf.y;
        if (IS3D) std::cout << " x " << blockConf.z;
        std::cout << "\ngrid:  " << gridConf.x << " x " << gridConf.y;
        if (IS3D) std::cout << " x " << gridConf.z;
        std::cout << std::endl;
    }
    
    for (int lb = 0; lb < input.lnblocks; lb++)
    {
        K_Init<<<gridConf, blockConf>>>(flow, err, input, lb);
        CuCheck(cudaPeekAtLastError());
    }
    CuCheck(cudaDeviceSynchronize());
}


#define stencilIdx(v,j) ((v)+(5+DIM)*(j))
#define f_DivSplit(q,j,l,v1)         (0.500*(q[stencilIdx((v1),(j))] + q[stencilIdx((v1),(j)+(l))]))
#define fg_QuadSplit(q,j,l,v1,v2)    (0.250*(q[stencilIdx((v1),(j))] + q[stencilIdx((v1),(j)+(l))])*(q[stencilIdx((v2),(j))] + q[stencilIdx((v2),(j)+(l))]))
#define fg_CubeSplit(q,j,l,v1,v2,v3) (0.125*(q[stencilIdx((v1),(j))] + q[stencilIdx((v1),(j)+(l))])*(q[stencilIdx((v2),(j))] + q[stencilIdx((v2),(j)+(l))])*(q[stencilIdx((v3),(j))] + q[stencilIdx((v3),(j)+(l))]))
#define fg_DivSplit(q,j,l,v1,v2)     (0.500*(q[stencilIdx((v1),(j)+(l))]*q[stencilIdx((v2),(j))]) + (q[stencilIdx((v1),(j))]*q[stencilIdx((v2),(j)+(l))]))

__global__ void K_Conv(double* flow, double* err, const InputClass input, const int lb, const Coef_t center, int stencilWid)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;
#if(IS3D)    
    int k = threadIdx.z + blockIdx.z*blockDim.z;
#else
    int k = 0;
#endif
    double invdx[DIM];
    for (int d = 0; d < DIM; d++) invdx[d] = input.nxb[d]/(input.bounds[2*d+1] - input.bounds[2*d]);
    int dijk[3] = {0};
    if (i < (input.nxb[0]) && j < (input.nxb[1]) && k < (input.nxb[2]))
    {
        for (int idir = 0; idir < DIM; idir++)
        {
            dijk[idir] = 1;
            double stencilData[9*(5+DIM)]; //ie,ke,rho,P,T,u,v,w
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
                double al = center.c[l-1];
                int jf = stencilWid;
                for (int m = 0; m <= (l-1); m++)
                {
                    C[1] += 2.0*al*fg_QuadSplit(stencilData,jf-m, l,2,5+idir);
                    C[0] += 2.0*al*fg_QuadSplit(stencilData,jf+m,-l,2,5+idir);

                    for (int idir_mom = 0; idir_mom < DIM; idir_mom++)
                    {
                        M[idir_mom      ] = 2.0*al*fg_CubeSplit(stencilData,jf-m, l,3,5+idir,5+idir_mom);
                        M[idir_mom + DIM] = 2.0*al*fg_CubeSplit(stencilData,jf+m,-l,3,5+idir,5+idir_mom);
                    }

                    PGRAD[1] += 2.0*al*f_DivSplit(stencilData,jf-m, l,3);
                    PGRAD[0] += 2.0*al*f_DivSplit(stencilData,jf+m,-l,3);

                    for (int vel_comp = 0;  vel_comp < DIM; vel_comp ++)
                    {
                        KE[1] += 2.0*al*fg_QuadSplit(stencilData,jf-m, l,2,5+idir)*0.5*(stencilData[stencilIdx(5+vel_comp,jf-m)]*stencilData[stencilIdx(5+vel_comp,jf-m+l)]);
                        KE[0] += 2.0*al*fg_QuadSplit(stencilData,jf+m,-l,2,5+idir)*0.5*(stencilData[stencilIdx(5+vel_comp,jf+m)]*stencilData[stencilIdx(5+vel_comp,jf+m-l)]);
                    }

                    IE[1] += 2.0*al*fg_CubeSplit(stencilData,jf-m, l,3,0,5+idir);
                    IE[0] += 2.0*al*fg_CubeSplit(stencilData,jf+m,-l,3,0,5+idir);

                    PDIFF[1] += 2.0*al*fg_DivSplit(stencilData,jf-m, l,5+idir,3);
                    PDIFF[0] += 2.0*al*fg_DivSplit(stencilData,jf-m, l,5+idir,3);
                }
            }
            double rhs[2+DIM] = {0.0};
            rhs[0] += -invdx[idir]*(C[1] - C[0]);
            rhs[1] += -invdx[idir]*(IE[1] + KE[1] + PDIFF[1] - IE[0] - KE[0] - PDIFF[0]);
            rhs[2+idir] += -invdx[idir]*(PGRAD[1] - PGRAD[0]);
            for (int rhs_vel_comp = 0; rhs_vel_comp < DIM; rhs_vel_comp++)
            {
                rhs[2+rhs_vel_comp] += -invdx[idir]*(M[1] - M[0]);
            }
            for (int myIndex = 0; myIndex < 2+DIM; myIndex++) err[bidx(0, i, j, k, lb, input)] = rhs[myIndex];
            dijk[idir] = 0;
        }
    }
}

void ConvGpu(double* flow, double* err, const InputClass& input)
{
    dim3 blockConf;
    blockConf.x = BLOCK_SIZE;
    blockConf.y = BLOCK_SIZE;
#if(IS3D)
    blockConf.z = BLOCK_SIZE;
#endif
    dim3 gridConf;
    int numcells[DIM];
    for (int i = 0; i < DIM; i++) {numcells[i] = input.nxb[i];}
    gridConf.x = (numcells[0] + BLOCK_SIZE - 1)/BLOCK_SIZE;
    gridConf.y = (numcells[1] + BLOCK_SIZE - 1)/BLOCK_SIZE;
#if(IS3D)
    gridConf.z = (numcells[2] + BLOCK_SIZE - 1)/BLOCK_SIZE;
#endif

    if (mypenoG==0 && !hasPrintedGp)
    {
        hasPrintedGp = true;
        std::cout << "GP Config:\nblock: " << blockConf.x << " x " << blockConf.y;
        if (IS3D) std::cout << " x " << blockConf.z;
        std::cout << "\ngrid:  " << gridConf.x << " x " << gridConf.y;
        if (IS3D) std::cout << " x " << gridConf.z;
        std::cout << std::endl;
    }
    
    Coef_t center;
    switch (input.centOrder)
    {
        case 2: {center.c[0] = 1.0/2.0; break;}
        case 4: {center.c[0] = 2.0/3.0; center.c[1] = -1.0/12.0; break;}
        case 6: {center.c[0] = 3.0/4.0; center.c[1] = -3.0/20.0; center.c[2] = 1.0/60.0 ;break;}
        case 8: {center.c[0] = 4.0/5.0; center.c[1] = -1.0/5.0 ; center.c[2] = 4.0/105; center.c[3] = -1.0/280.0; break;}
        default: {std::cout << "Bad central scheme order." << std::endl; abort();}
    }
    int stencilWid = input.centOrder/2;
    for (int lb = 0; lb < input.lnblocks; lb++)
    {
        K_Conv<<<gridConf, blockConf>>>(flow, err, input, lb, center, stencilWid);
        CuCheck(cudaPeekAtLastError());
    }
    CuCheck(cudaDeviceSynchronize());


}