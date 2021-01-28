#ifndef IDX_H
#define IDX_H

#ifdef __CUDACC__
#define GPU 1
#else
#define GPU 0
#endif

#if (GPU)



#else



#endif

#endif