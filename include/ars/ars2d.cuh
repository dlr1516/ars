#include <ars/functions.h>

//#include <ars/ars2d.h>

// --------------------------------------------------------
// 2D->1D INDICIZATION IN FOURIER COEFFICIENT MATRIX
// --------------------------------------------------------

/**
 * When dealing with Fourier coefficient matrix, return index referring to the first point that is being dealt with.
 * In short, when computing Ars(means[i], means[j]) -> this function returns i.
 */
__device__
int getIfromTid(int tid, int n);

/**
 * When dealing with Fourier coefficient matrix, return index referring to the second point that is being dealt with.
 * In short, when computing Ars(means[i], means[j]) -> this function returns j.
 */
__device__
int getJfromTid(int tid, int n, int i);

// --------------------------------------------------------
// PNEBI FUNCTIONS
// PNEBI stands for Product of Negative Exponential and Bessel I, which is defined as
// 
//    PNEBI(k,x) = 2.0 * exp(-x) * besseli(k,x)
// 
// where besseli(k,x) is the modified Bessel function of the First Kind with order k. 
// --------------------------------------------------------

/** 
 * Computes the value of function of PNEBI(0,x) = 2.0 * exp(-x) * besseli(0,x) 
 * using the polynomial approximation of Abramowitz-Stegun (9.8.1)-(9.8.2). 
 * Common library functions computing besseli(0,x) leads to numeric overflow 
 * or exploits inaccurate (and substantially flat in our interval!) Hankel 
 * approximation. 
 */
__device__
double evaluatePnebi0Polynom(double x);

/** 
 * Evaluates PNEBI function in point x for different orders from 0 to n. 
 * This implementation is based on downward recurring formula as suggested in
 *  
 * Aa Vv, Numerical Recipes in C. The Art of Scientific Computing, 2nd edition, 1992. 
 * 
 * It uses down
 */
__device__
void evaluatePnebiVectorGPU(int n, double x, double* pnebis, int pnebisSz);

// --------------------------------------------------------
// GLOBAL CUDA KERNELS
// --------------------------------------------------------

/**
 * Insert Isotropic Gaussians Kernel, that uses Downward method for partial coefficients computing
 * @param means
 * @param sigma1
 * @param sigma2
 * @param numPts
 * @param fourierOrder
 * @param numColsPadded
 * @param pnebiMode
 * @param coeffsMat
 */
__global__
void iigDw(cuars::Vec2d* means, double sigma1, double sigma2, int numPts, int fourierOrder, int numColsPadded, cuars::ArsKernelIso2dComputeMode pnebiMode, double* coeffsMat);

/**
 * !! UNFINISHED
 * Insert Isotropic Gaussians Kernel, that uses LUT table method for partial coefficients computing
 * @param means
 * @param sigma1
 * @param sigma2
 * @param numPts
 * @param numPtsAfterPadding
 * @param fourierOrder
 * @param numColsPadded
 * @param pnebiMode
 * @param pnebiLUT
 * @param coeffsMat
 */
__global__
void iigLut(cuars::Vec2d* means, double sigma1, double sigma2, int numPts, int numPtsAfterPadding, int fourierOrder, int numColsPadded, cuars::ArsKernelIso2dComputeMode pnebiMode, cuars::PnebiLUT& pnebiLUT, double* coeffsMat);

__global__
void makePartialSums(double* matIn, int nrowsIn, int ncols, double *matOut);

__global__
void sumColumnsPartialSums(double* matIn, int nrows, int ncols, double* vecOut);