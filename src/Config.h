#ifndef CONFIG_H
#define CONFIG_H

#ifndef DIM
#define DIM 3
#endif

#if (DIM == 3)
#define IS3D 1
#else
#define IS3D 0
#endif

#define BLOCK_SIZE 8

#endif