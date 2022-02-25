#include <ars/definitions.h>
#include <ars/ars2d.h>
#include <ars/BBOptimizer1d.h>

#include <thrust/host_vector.h>
#include <thrust/device_malloc.h>

int ceilPow2(int n);

int sumNaturalsUpToN(int n);

__device__
int getIfromTid(int tid, int n);

__device__
int getJfromTid(int tid, int n, int i);

__device__
double evaluatePnebi0Polynom(double x);

__device__
void evaluatePnebiVectorGPU(int n, double x, double* pnebis, int pnebisSz);

__global__
void iigKernelDownward_old(cuars::Vec2d* means, double sigma1, double sigma2, int numPts, int numPtsAfterPadding, int fourierOrder, int numColsPadded, cuars::ArsKernelIsotropic2d::ComputeMode pnebiMode, double* coeffsMat);

__global__
void iigKernelDownward(cuars::Vec2d* means, double sigma1, double sigma2, int numPts, int fourierOrder, int numColsPadded, cuars::ArsKernelIsotropic2d::ComputeMode pnebiMode, double* coeffsMat);

__global__
void iigKernelLut(cuars::Vec2d* means, double sigma1, double sigma2, int numPts, int numPtsAfterPadding, int fourierOrder, int numColsPadded, cuars::ArsKernelIsotropic2d::ComputeMode pnebiMode, cuars::PnebiLUT& pnebiLUT, double* coeffsMat);

__global__
void sumColumns(double* mat, int nrows, int ncols, double* sums);

__global__
void sumColumnsNoPadding(double* mat, int nrows, int ncols, double* sums);