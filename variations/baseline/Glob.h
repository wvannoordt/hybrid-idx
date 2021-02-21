#ifndef GLOB_H
#define GLOB_H

#include <string>

extern int mypenoG;
extern bool hasPrintedGp;

inline std::string zfill(int n, int numSpaces)
{
    std::string output = std::to_string(n);
    while (output.length() < numSpaces) output = "0" + output;
    return output;
}

#endif