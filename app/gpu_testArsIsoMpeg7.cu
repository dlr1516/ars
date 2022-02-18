
/**
 * ARS - Angular Radon Spectrum 
 * Copyright (C) 2017 Dario Lodi Rizzini.
 *
 * ARS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * ARS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with ARS.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <iostream>
#include <chrono>
#include <cuda_runtime.h>

#include "ars/mpeg7RW.h"
#include "ars/cuArsIso.cuh"



#define PRINT_DIM(X) std::cout << #X << " rows " << X.rows() << " cols " << X.cols() << std::endl;
#define RAD2DEG(X) (180.0/M_PI*(X))

struct BoundInterval {
    double x0;
    double x1;
    double y0;
    double y1;
};

int main(void) {

    ArsImgTests::PointReaderWriter pointsSrc;
    ArsImgTests::PointReaderWriter pointsDst;

    std::string filenameSrc = "/home/rimlab/Downloads/mpeg7_point_tests/noise000_occl00_rand000/apple-1_xp0686_yp0967_t059_sigma0001_occl000.txt";
    std::string filenameDst = "/home/rimlab/Downloads/mpeg7_point_tests/noise000_occl00_rand000/apple-1_xp0749_yn0521_t090_sigma0001_occl000.txt";

    // Loads files and computes the rotation
    std::cout << "\n*****\nLoading file \"" << filenameSrc << "\"" << std::endl;
    pointsSrc.load(filenameSrc);
    std::cout << "\n*****\nLoading file \"" << filenameDst << "\"" << std::endl;
    pointsDst.load(filenameDst);
    std::cout << "  points src " << pointsSrc.points().size() << ", points dst " << pointsDst.points().size() << std::endl;



    cuars::AngularRadonSpectrum2d arsSrc;
    cuars::AngularRadonSpectrum2d arsDst;
    double sigma = 0.05;
    int fourierOrder = 20;

    //ARS parameters setting
    arsSrc.setARSFOrder(fourierOrder);
    arsDst.setARSFOrder(fourierOrder);
    cuars::ArsKernelIsotropic2d::ComputeMode pnebiMode = cuars::ArsKernelIsotropic2d::ComputeMode::PNEBI_DOWNWARD;
    arsSrc.setComputeMode(pnebiMode);
    arsDst.setComputeMode(pnebiMode);


    //Fourier coefficients mega-matrix computation -> parallelization parameters
    int numPts = std::min<int>(pointsSrc.points().size(), pointsDst.points().size()); //the two should normally be equal
    const int numPtsAfterPadding = ceilPow2(numPts); //for apple1 -> numPts 661; padded 1024
    const int blockSize = 256; //num threads per block
    const int numBlocks = (numPtsAfterPadding * numPtsAfterPadding) / blockSize; //number of blocks in grid (each block contains blockSize threads)
    const int gridTotalSize = blockSize*numBlocks; //total number of threads in grid
    //depth of mega-matrix
    const int coeffsMatNumCols = 2 * fourierOrder + 2;
    const int coeffsMatNumColsPadded = ceilPow2(coeffsMatNumCols);
    const int coeffsMatTotalSz = numPtsAfterPadding * numPtsAfterPadding * coeffsMatNumColsPadded;
    //Fourier matrix sum -> parallelization parameters
    const int sumBlockSz = 64;
    const int sumGridSz = 256; //can be used to futher parallelize sum of mega-matrix (for now in sum kernel it is actually set to 1)
    std::cout << "Parallelization params:" << std::endl;
    std::cout << "numPtsAfterPadding " << numPtsAfterPadding << " blockSize " << blockSize << " numBlocks " << numBlocks << " gridTotalSize " << gridTotalSize << std::endl;
    std::cout << "sumBlockSz " << sumBlockSz << " sumGridSz " << sumGridSz << std::endl;

    std::cout << "\n------\n" << std::endl;

    std::cout << "\n\nCalling kernel functions on GPU\n" << std::endl;


    //ARS SRC -> preparation for kernel calls and kernel calls
    cudaEvent_t startSrc, stopSrc; //timing using CUDA events
    cudaEventCreate(&startSrc);
    cudaEventCreate(&stopSrc);
    cuars::Vec2d * kernelInputSrc;
    cudaMalloc((void**) &kernelInputSrc, numPtsAfterPadding * sizeof (cuars::Vec2d));
    cudaMemcpy(kernelInputSrc, pointsSrc.points().data(), numPtsAfterPadding * sizeof (cuars::Vec2d), cudaMemcpyHostToDevice);

    double *coeffsMatSrc;
    cudaMalloc((void**) &coeffsMatSrc, coeffsMatTotalSz * sizeof (double));
    cudaMemset(coeffsMatSrc, 0.0, coeffsMatTotalSz * sizeof (double));
    //    for (int i = 0; i < coeffsMatTotalSz; ++i) {
    //        coeffsMaSrc1[i] = 0.0;
    //    }
    double* d_coeffsArsSrc;
    cudaMalloc((void**) &d_coeffsArsSrc, coeffsMatNumColsPadded * sizeof (double));
    cudaMemset(d_coeffsArsSrc, 0.0, coeffsMatNumColsPadded * sizeof (double));

    cudaEventRecord(startSrc);
    iigKernelDownward << <numBlocks, blockSize >> >(kernelInputSrc, sigma, sigma, numPts, numPtsAfterPadding, fourierOrder, coeffsMatNumColsPadded, pnebiMode, coeffsMatSrc);
    sumColumns << <1, sumBlockSz>> >(coeffsMatSrc, numPtsAfterPadding, coeffsMatNumColsPadded, d_coeffsArsSrc);
    cudaEventRecord(stopSrc);

    double* coeffsArsSrc = new double [coeffsMatNumColsPadded];
    cudaMemcpy(coeffsArsSrc, d_coeffsArsSrc, coeffsMatNumColsPadded * sizeof (double), cudaMemcpyDeviceToHost);

    cudaEventSynchronize(stopSrc);
    float millisecondsSrc = 0.0f;
    cudaEventElapsedTime(&millisecondsSrc, startSrc, stopSrc);
    std::cout << "SRC -> insertIsotropicGaussians() " << millisecondsSrc << " ms" << std::endl;

    cudaError_t cudaerr = cudaDeviceSynchronize();
    if (cudaerr != cudaSuccess)
        printf("kernel launch failed with error \"%s\".\n", cudaGetErrorString(cudaerr));
    
    //    for (int i = 0; i < coeffsMatNumColsPadded; ++i) {
    //        std::cout << "coeffsArsSrc[" << i << "] " << coeffsArsSrc[i] << std::endl;
    //    }

    cudaFree(coeffsMatSrc);
    cudaFree(kernelInputSrc);
    cudaFree(d_coeffsArsSrc);
    cudaEventDestroy(startSrc);
    cudaEventDestroy(stopSrc);
    //END OF ARS SRC



    std::cout << "\n------\n" << std::endl; //"pause" between ars src and ars dst



    //ARS DST -> preparation for kernel calls and kernel calls
    cudaEvent_t startDst, stopDst; //timing using CUDA events
    cudaEventCreate(&startDst);
    cudaEventCreate(&stopDst);
    cuars::Vec2d *kernelInputDst;
    cudaMalloc((void**) &kernelInputDst, numPtsAfterPadding * sizeof (cuars::Vec2d));
    cudaMemcpy(kernelInputDst, pointsDst.points().data(), numPtsAfterPadding * sizeof (cuars::Vec2d), cudaMemcpyHostToDevice);

    double *coeffsMatDst; //magari evitare di fare il delete e poi riallocarla è più efficiente (anche se comunque ci sarebbe poi da settare tutto a 0)
    cudaMalloc((void**) &coeffsMatDst, coeffsMatTotalSz * sizeof (double));
    cudaMemset(coeffsMatDst, 0.0, coeffsMatTotalSz * sizeof (double));
    //    for (int i = 0; i < coeffsMatTotalSz; ++i) {
    //        coeffsMatDst[i] = 0.0;
    //    }
    double* d_coeffsArsDst;
    cudaMalloc((void**) &d_coeffsArsDst, coeffsMatNumColsPadded * sizeof (double));
    cudaMemset(d_coeffsArsDst, 0.0, coeffsMatNumColsPadded * sizeof (double));

    cudaEventRecord(startDst);
    iigKernelDownward << <numBlocks, blockSize >> >(kernelInputDst, sigma, sigma, numPts, numPtsAfterPadding, fourierOrder, coeffsMatNumColsPadded, pnebiMode, coeffsMatDst);
    sumColumns << <1, sumBlockSz>> >(coeffsMatDst, numPtsAfterPadding, coeffsMatNumColsPadded, d_coeffsArsDst);
    cudaEventRecord(stopDst);



    double* coeffsArsDst = new double [coeffsMatNumColsPadded];
    cudaMemcpy(coeffsArsDst, d_coeffsArsDst, coeffsMatNumColsPadded * sizeof (double), cudaMemcpyDeviceToHost);

    cudaEventSynchronize(stopDst);
    float millisecondsDst = 0.0f;
    cudaEventElapsedTime(&millisecondsDst, startDst, stopDst);
    std::cout << "DST -> insertIsotropicGaussiansDst() " << millisecondsDst << " ms" << std::endl;

    cudaerr = cudaDeviceSynchronize();
    if (cudaerr != cudaSuccess)
        printf("kernel launch failed with error \"%s\".\n", cudaGetErrorString(cudaerr));
    
    //    for (int i = 0; i < coeffsMatNumColsPadded; ++i) {
    //        std::cout << "coeffsArsDst[" << i << "] " << coeffsArsDst[i] << std::endl;
    //    }

    cudaFree(coeffsMatDst);
    cudaFree(kernelInputDst);
    cudaFree(d_coeffsArsDst);
    cudaEventDestroy(startDst);
    cudaEventDestroy(stopDst);
    //END OF ARS DST



    //Computation final computations (correlation, ...) on CPU
    std::cout << "\nARS Coefficients:\n";
    std::cout << "\ti \tDownward \tLUT\n";
    arsSrc.setCoefficients(coeffsArsSrc, coeffsMatNumCols);
    //    for (int i = 0; i < coeffsVectorMaxSz; i++) {
    //        std::cout << "arsSrc - coeff_d[" << i << "] " << d_coeffsMat1[i] << std::endl;
    //    }
    arsDst.setCoefficients(coeffsArsDst, coeffsMatNumCols);
    for (int i = 0; i < arsSrc.coefficients().size() && i < arsDst.coefficients().size(); ++i) {
        std::cout << "\t" << i << " \t" << arsSrc.coefficients().at(i) << " \t" << arsDst.coefficients().at(i) << "\n";
    }
    std::cout << std::endl;

    std::vector<double> funcFourierRecursDownSrc;
    std::vector<double> funcFourierRecursDownDst;
    int thnum = 360;
    double dtheta = M_PI / thnum;
    double theta;
    for (int i = 0; i < thnum; ++i) {
        theta = dtheta * i;
        funcFourierRecursDownSrc.push_back(arsSrc.eval(theta));
        funcFourierRecursDownDst.push_back(arsDst.eval(theta));
    }

    std::cout << "\nBranch and Bound limits:\n";
    int bbnum = 32;
    std::vector<BoundInterval> bbbs(bbnum);
    for (int i = 0; i < bbnum; ++i) {
        bbbs[i].x0 = M_PI * i / bbnum;
        bbbs[i].x1 = M_PI * (i + 1) / bbnum;
        cuars::findLUFourier(arsSrc.coefficients(), bbbs[i].x0, bbbs[i].x1, bbbs[i].y0, bbbs[i].y1);
        std::cout << i << ": x0 " << RAD2DEG(bbbs[i].x0) << " x1 " << RAD2DEG(bbbs[i].x1) << ", y0 " << bbbs[i].y0 << " y1 " << bbbs[i].y1 << std::endl;
    }


    cuars::FourierOptimizerBB1D optim(arsSrc.coefficients());
    double xopt, ymin, ymax;
    optim.enableXTolerance(true);
    optim.enableYTolerance(true);
    optim.setXTolerance(M_PI / 180.0 * 0.5);
    optim.setYTolerance(1.0);
    optim.findGlobalMax(0, M_PI, xopt, ymin, ymax);
    std::cout << "\n****\nMaximum in x = " << xopt << " (" << RAD2DEG(xopt) << " deg), maximum between [" << ymin << "," << ymax << "]" << std::endl;

    double xopt2, ymax2;
    cuars::findGlobalMaxBBFourier(arsSrc.coefficients(), 0, M_PI, M_PI / 180.0 * 0.5, 1.0, xopt2, ymax2);
    std::cout << "  repeated evaluation with findGlobalMaxBBFourier(): maximum in x " << xopt2 << " (" << RAD2DEG(xopt2) << " deg), maximum value " << ymax2 << std::endl;



    //Free CPU memory
    free(coeffsArsSrc);
    free(coeffsArsDst);


    return 0;
}





