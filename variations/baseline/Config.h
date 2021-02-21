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

#define BLOCK_SIZEX 4
#define BLOCK_SIZEY 4
#define BLOCK_SIZEZ 4

#endif