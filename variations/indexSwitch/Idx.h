#ifndef IDX_H
#define IDX_H
#include "Config.h"
#ifdef __CUDACC__
#define GPU 1
#else
#define GPU 0
#endif

#if (GPU)
// ~70 ms
#define bidx(var,i,j,k,lb,input) \
(\
         ((i)+input.nguard) + \
         ((j)+input.nguard)*(input.nxb[0]+2*input.nguard) + \
    IS3D*((k)+input.nguard)*(input.nxb[0]+2*input.nguard)*(input.nxb[1]+2*input.nguard) + \
         (              lb)*(input.nxb[0]+2*input.nguard)*(input.nxb[1]+2*input.nguard)*(IS3D*input.nxb[1+IS3D]+2*input.nguard + (1-IS3D)) + \
         (             var)*(input.nxb[0]+2*input.nguard)*(input.nxb[1]+2*input.nguard)*(IS3D*input.nxb[1+IS3D]+2*input.nguard + (1-IS3D)) * (input.lnblocks) \
)

// ~70 ms
// #define bidx(var,i,j,k,lb,input) \
// (\
//     (i) + input.nguard + \
//     ((j)+input.nguard)*(input.nxb[0]+2*input.nguard) + \
//     IS3D*((k)+input.nguard)*(input.nxb[0]+2*input.nguard)*(input.nxb[1]+2*input.nguard) + \
//          (              lb)*(input.nxb[0]+2*input.nguard)*(input.nxb[1]+2*input.nguard)*(IS3D*input.nxb[1+IS3D]+2*input.nguard + (1-IS3D)) + \
//     (var) * (input.lnblocks)*(input.nxb[0]+2*input.nguard)*(input.nxb[1]+2*input.nguard)*(IS3D*input.nxb[1+IS3D]+2*input.nguard + (1-IS3D))\
// )

#else

// ~122 ms
#define bidx(var,i,j,k,lb,input) \
(\
     ((i)+input.nguard) + \
     ((j)+input.nguard)*(input.nxb[0]+2*input.nguard) + \
IS3D*((k)+input.nguard)*(input.nxb[0]+2*input.nguard)*(input.nxb[1]+2*input.nguard) + \
     (              lb)*(input.nxb[0]+2*input.nguard)*(input.nxb[1]+2*input.nguard)*(IS3D*input.nxb[1+IS3D]+2*input.nguard + (1-IS3D)) + \
     (             var)*(input.nxb[0]+2*input.nguard)*(input.nxb[1]+2*input.nguard)*(IS3D*input.nxb[1+IS3D]+2*input.nguard + (1-IS3D)) * (input.lnblocks) \
)

// 122~ ms
// #define bidx(var,i,j,k,lb,input) \
// (\
//     (i) + input.nguard + \
//     ((j)+input.nguard)*(input.nxb[0]+2*input.nguard) + \
//     IS3D*((k)+input.nguard)*(input.nxb[0]+2*input.nguard)*(input.nxb[1]+2*input.nguard) + \
//          (              lb)*(input.nxb[0]+2*input.nguard)*(input.nxb[1]+2*input.nguard)*(IS3D*input.nxb[1+IS3D]+2*input.nguard + (1-IS3D)) + \
//     (var) * (input.lnblocks)*(input.nxb[0]+2*input.nguard)*(input.nxb[1]+2*input.nguard)*(IS3D*input.nxb[1+IS3D]+2*input.nguard + (1-IS3D))\
// )

#endif

#endif