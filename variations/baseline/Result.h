#ifndef HI_RESULT_H
#define HI_RESULT_H
#include <fstream>
#include <string>
struct Result
{
    double avgStepTime;
    void WriteJson(std::ofstream& fs)
    {
        fs << "    \"avgStepTime\": " << avgStepTime << "\n";
    }
};

#endif