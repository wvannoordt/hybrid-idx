#ifndef IDX_H
#define IDX_H
#include "Config.h"
#ifdef __CUDACC__
#define GPU 1
#else
#define GPU 0
#endif

#if (GPU)



#else

#define bidx(var,i,j,k,lb,input) \
(\
    var + \
         ((i)+input.nguard)*(2+DIM) + \
         ((j)+input.nguard)*(input.nxb[0]+2*input.nguard)*(2+DIM) + \
    IS3D*((k)+input.nguard)*(input.nxb[0]+2*input.nguard)*(input.nxb[1]+2*input.nguard)*(2+DIM) + \
         (              lb)*(input.nxb[0]+2*input.nguard)*(input.nxb[1]+2*input.nguard)*(IS3D*input.nxb[1+IS3D]+2*input.nguard + (1-IS3D))*(2+DIM)\
)

#endif

#endif