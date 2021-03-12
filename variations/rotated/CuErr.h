#ifndef CU_ERR_H
#define CU_ERR_H
#include <iostream>
#define CuCheck(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      std::cout << "Failed CUDA call in file " << file << ", line " << line << ". Message:\n" << cudaGetErrorString(code) << std::endl;
      if (abort) exit(code);
   }
}

#endif