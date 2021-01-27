#include "PTL.h"
#include "mpi.h"
#include "Config.h"
struct Inputclass
{
    int centOrder;
    int* nxb;
    int lnblocks;
} input;

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int mypeno, noprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &mypeno);
    MPI_Comm_size(MPI_COMM_WORLD, &noprocs);
    PTL::PropertyTree ptree;
    ptree["centOrder"].MapTo(&input.centOrder) = new PTL::PTLInteger(2, "order of central scheme");
    ptree["nxb"].MapTo(&input.nxb) = new PTL::PTLStaticIntegerArray(DIM, "dimensions of blocks", [](int i){return 16;});
    ptree["lnblocks"].MapTo(&input.lnblocks) = new PTL::PTLInteger(10, "total number of blocks");
    ptree.Read("input.ptl");
    ptree.StrictParse();
    
    std::cout << "Hello from " << mypeno << " of " << noprocs << std::endl;
    
    MPI_Finalize();
    
}