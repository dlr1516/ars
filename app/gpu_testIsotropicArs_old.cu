
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

#include "ars/cuArsIso.cuh"



#define PRINT_DIM(X) std::cout << #X << " rows " << X.rows() << " cols " << X.cols() << std::endl;
#define RAD2DEG(X) (180.0/M_PI*(X))

struct BoundInterval {
    double x0;
    double x1;
    double y0;
    double y1;
};

void rangeToPoint(double* ranges, int num, int numPadded, double angleMin, double angleRes, std::vector<cuars::Vec2d>& points);

int main(void) {
    double acesRanges[] = {50.00, 50.00, 50.00, 5.26, 5.21, 5.06, 5.01, 3.01, 2.94, 2.89, 2.84, 2.74, 2.69, 2.64, 2.59, 2.54, 2.49, 2.49, 2.44, 2.39, 2.34, 2.29, 2.29, 2.29, 2.39, 2.39, 2.49, 2.51, 2.61, 2.66, 2.76, 2.81, 2.96, 3.01, 3.11, 3.26, 3.01, 3.01, 3.01, 3.06, 3.21, 6.86, 6.86, 6.81, 6.76, 6.71, 6.71, 6.66, 6.61, 6.66, 6.56, 6.56, 6.56, 6.46, 6.46, 6.41, 6.46, 6.46, 4.11, 3.96, 3.96, 4.96, 4.86, 5.21, 7.41, 4.61, 5.16, 6.26, 6.26, 6.31, 4.86, 5.01, 5.86, 5.81, 4.21, 4.26, 4.31, 4.41, 4.39, 4.46, 5.31, 5.06, 5.26, 4.96, 6.01, 5.76, 5.61, 5.36, 5.26, 5.01, 4.21, 4.16, 4.01, 3.91, 3.61, 3.21, 3.26, 3.16, 3.06, 3.01, 3.31, 3.21, 3.16, 2.16, 2.19, 2.16, 2.21, 2.11, 2.01, 2.01, 2.06, 2.84, 2.91, 2.91, 3.01, 3.11, 3.21, 3.81, 4.06, 7.11, 7.06, 7.01, 6.96, 6.86, 4.31, 6.76, 6.71, 6.66, 6.61, 5.46, 5.41, 6.46, 6.21, 6.31, 6.51, 7.26, 7.46, 50.00, 2.01, 1.94, 1.94, 1.94, 2.31, 1.86, 1.84, 1.84, 1.81, 1.96, 26.46, 20.76, 2.11, 2.12, 2.17, 2.14, 2.09, 2.09, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.19, 2.19, 2.24, 2.24, 2.24, 2.24, 2.29, 2.29, 2.29, 2.29, 2.29, 2.39, 2.39, 2.39, 2.44};
    cuars::AngularRadonSpectrum2d ars1;
    cuars::AngularRadonSpectrum2d ars2;
    //    std::chrono::system_clock::time_point timeStart, timeStop;
    double sigma = 0.05;
    int fourierOrder = 20;

    ars1.setARSFOrder(fourierOrder);
    ars2.setARSFOrder(fourierOrder);

    //parallelization parameters
    int numPts = 180; // = acesRanges.size()
    const int numPtsAfterPadding = ceilPow2(numPts);
    const int blockSize = 256; //num threads per block
    const int numBlocks = (numPtsAfterPadding * numPtsAfterPadding) / blockSize; //number of blocks in grid (each block contains blockSize threads)
    const int gridTotalSize = blockSize*numBlocks; //total number of threads in grid

    const int sumBlockSz = 64;
    const int sumGridSz = 256;
    std::cout << "Parallelization params:" << std::endl;
    std::cout << "numPtsAfterPadding " << numPtsAfterPadding << " blockSize " << blockSize << " numBlocks " << numBlocks << " gridTotalSize " << gridTotalSize << std::endl;
    std::cout << "sumSrcBlockSz " << sumBlockSz << " sumGridSz " << sumGridSz << std::endl;
    std::cout << "numPtsAfterPadding " << numPtsAfterPadding << " blockSize " << blockSize << " numBlocks " << numBlocks << " gridTotalSize " << gridTotalSize << std::endl;

    //conversion
    std::vector<cuars::Vec2d> acesPointsSTL;
    //    cuars::Vec2d p0, p1;
    //    p0.x = 0.0;
    //    p0.y = 0.0;
    //    p1.x = cos(M_PI * 30 / 180.0);
    //    p1.y = sin(M_PI * 30 / 180.0);
    //    acesPointsSTL.push_back(p0);
    //    acesPointsSTL.push_back(p1);

    rangeToPoint(acesRanges, numPts, numPtsAfterPadding, -0.5 * M_PI, M_PI / 180.0 * 1.0, acesPointsSTL);

    thrust::host_vector<cuars::Vec2d> acesPointsHost(acesPointsSTL.begin(), acesPointsSTL.end());


    //    cuars::Vec2d firstElement; //??

    std::cout << "\n------\n" << std::endl;
    std::cout << "\n\nCalling kernel functions on GPU\n" << std::endl;

    //    timeStart = std::chrono::system_clock::now();
    //ars1 kernel call
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    //    ars1.insertIsotropicGaussians(acesPoints1, sigma);
    cuars::Vec2d * kernelInput1;
    cudaMalloc((void**) &kernelInput1, numPtsAfterPadding * sizeof (cuars::Vec2d));
    cudaMemcpy(kernelInput1, acesPointsHost.data(), numPtsAfterPadding * sizeof (cuars::Vec2d), cudaMemcpyHostToDevice);
    //    for (int i = 0; i < numPtsAfterPadding; ++i) {
    //        kernelInput1[i] = acesPointsSTL[i];
    //    }

    //    cudaDeviceSynchronize();
    //    std::cout << "acesPointsHost.size() " << acesPointsHost.size() << std::endl;
    //    for (int s = 0; s < acesPointsHost.size(); s++) {
    //        std::cout << "s " << s << std::endl;
    //        std::cout << kernelInput1[s].x << " " << kernelInput1[s].y << std::endl;
    //    }


    //    ars1.initLUT(0.0001);
    cuars::ArsKernelIsotropic2d::ComputeMode pnebiMode = cuars::ArsKernelIsotropic2d::ComputeMode::PNEBI_DOWNWARD;
    ars1.setComputeMode(pnebiMode);


    const int coeffsMatNumCols = 2 * fourierOrder + 2;
    const int coeffsMatNumColsPadded = ceilPow2(coeffsMatNumCols);
    const int coeffsMatTotalSz = numPtsAfterPadding * numPtsAfterPadding * coeffsMatNumColsPadded;
    double *coeffsMat1;
    cudaMalloc((void**) &coeffsMat1, coeffsMatTotalSz * sizeof (double));
    cudaMemset(coeffsMat1, 0.0, coeffsMatTotalSz * sizeof (double));
    //    for (int i = 0; i < coeffsMatTotalSz; ++i) {
    //        coeffsMat1[i] = 0.0;
    //    }
    double* d_coeffsArs1;
    cudaMalloc((void**) &d_coeffsArs1, coeffsMatNumColsPadded * sizeof (double));
    cudaMemset(d_coeffsArs1, 0.0, coeffsMatNumColsPadded * sizeof (double));


    //    cuars::PnebiLUT pnebiLUT1; //LUT setup
    //    double lutPrecision = 0.001; //LUT setup
    //    pnebiLUT1.init(fourierOrder, lutPrecision); //LUT setup
    //    if (pnebiLUT1.getOrderMax() < fourierOrder) { //LUT setup
    //        ARS_ERROR("LUT not initialized to right order. Initialized now."); //LUT setup
    //        pnebiLUT1.init(fourierOrder, 0.0001); //LUT setup
    //    }

    cudaEventRecord(start);
    iigKernelDownward_old << <numBlocks, blockSize >> >(kernelInput1, sigma, sigma, numPts, numPtsAfterPadding, fourierOrder, coeffsMatNumColsPadded, pnebiMode, coeffsMat1);
    sumColumns << <1, sumBlockSz>> >(coeffsMat1, numPtsAfterPadding, coeffsMatNumColsPadded, d_coeffsArs1);
    cudaEventRecord(stop);

    double* coeffsArs1 = new double [coeffsMatNumColsPadded];
    cudaMemcpy(coeffsArs1, d_coeffsArs1, coeffsMatNumColsPadded * sizeof (double), cudaMemcpyDeviceToHost);

    cudaEventSynchronize(stop);
    float millisecondsArs1 = 0.0;
    cudaEventElapsedTime(&millisecondsArs1, start, stop);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    //END OF ARS1

    std::cout << "\n------\n" << std::endl;


    //ARS2    
    ars2.setComputeMode(cuars::ArsKernelIsotropic2d::ComputeMode::PNEBI_LUT);

    cuars::PnebiLUT pnebiLUT2; //LUT setup
    double lutPrecision = 0.001; //already initialized for pnebiLUT1
    pnebiLUT2.init(fourierOrder, lutPrecision); //LUT setup
    if (pnebiLUT2.getOrderMax() < fourierOrder) { //LUT setup
        ARS_ERROR("LUT not initialized to right order. Initialized now."); //LUT setup
        pnebiLUT2.init(fourierOrder, 0.0001); //LUT setup
    }


    //kernel call
    //    timeStart = std::chrono::system_clock::now();
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    //    ars2.insertIsotropicGaussians(acesPoints1, sigma);

    cuars::Vec2d* kernelInput2;
    cudaMalloc((void **) &kernelInput2, numPtsAfterPadding * sizeof (cuars::Vec2d));
    pnebiMode = cuars::ArsKernelIsotropic2d::ComputeMode::PNEBI_LUT;

    double *coefficientsArs2 = new double[coeffsMatTotalSz](); //() initialize to 0
    double *d_coefficientsArs2; //d_ stands for device
    //    const int coeffsVectorMaxSz = 2 * fourierOrder + 2; //already initialized in ars1
    cudaMalloc(&d_coefficientsArs2, coeffsMatTotalSz * sizeof (double)); //maybe directly use cudaMemset?
    cudaMemcpy(d_coefficientsArs2, coefficientsArs2, coeffsMatTotalSz * sizeof (double), cudaMemcpyHostToDevice);


    cudaEventRecord(start);
    //    iigKernel << < numBlocks, blockSize >> >(thrust::raw_pointer_cast<ars::Vec2d*>(kernelInput2.data()), sigma, sigma, numPts, paddedPtVecSz, fourierOrder, pnebiMode, pnebiLUT2, d_coefficientsArs2);
    //    sumColumns << <1, sumBlockSz>> >(coeffsMat2, numPtsAfterPadding, coeffsMatNumColsPadded, d_coeffsArs2);
    cudaEventRecord(stop);
    //end of kernel calls for ARS2

    cudaMemcpy(coefficientsArs2, d_coefficientsArs2, coeffsMatTotalSz * sizeof (double), cudaMemcpyDeviceToHost);

    cudaEventSynchronize(stop);
    float millisecondsArs2 = 0.0;
    cudaEventElapsedTime(&millisecondsArs2, start, stop);

    //    timeStop = std::chrono::system_clock::now();
    //    double timeArs2 = (double) std::chrono::duration_cast<std::chrono::milliseconds>(timeStop - timeStart).count();
    //    std::cout << "insertIsotropicGaussians() " << timeArs2 << " ms" << std::endl;
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    //END OF ARS2


    cudaError_t cudaerr = cudaDeviceSynchronize();
    if (cudaerr != cudaSuccess)
        printf("kernel launch failed with error \"%s\".\n",
            cudaGetErrorString(cudaerr));

    std::cout << "ARS1 execution " << millisecondsArs1 << " ms --- ARS2 execution " << millisecondsArs2 << " ms" << std::endl;



    std::cout << "\nARS Coefficients:\n";
    std::cout << "\ti \tDownward \tLUT\n";
    ars1.setCoefficients(coeffsArs1, coeffsMatNumCols);
    //    for (int i = 0; i < coeffsVectorMaxSz; i++) {
    //        std::cout << "ars1coeff_d[" << i << "] " << d_coeffsMat1[i] << std::endl;
    //    }
    ars2.setCoefficients(coefficientsArs2, coeffsMatNumCols);
    for (int i = 0; i < ars1.coefficients().size() && i < ars2.coefficients().size(); ++i) {
        std::cout << "\t" << i << " \t" << ars1.coefficients().at(i) << " \t" << ars2.coefficients().at(i) << "\n";
    }
    std::cout << std::endl;

    std::vector<double> funcFourierRecursDownLUT;
    std::vector<double> funcFourierRecursDown;
    int thnum = 360;
    double dtheta = M_PI / thnum;
    double theta;
    for (int i = 0; i < thnum; ++i) {
        theta = dtheta * i;
        funcFourierRecursDownLUT.push_back(ars1.eval(theta));
        funcFourierRecursDown.push_back(ars2.eval(theta));
    }

    std::cout << "\nBranch and Bound limits:\n";
    int bbnum = 32;
    std::vector<BoundInterval> bbbs(bbnum);
    for (int i = 0; i < bbnum; ++i) {
        bbbs[i].x0 = M_PI * i / bbnum;
        bbbs[i].x1 = M_PI * (i + 1) / bbnum;
        cuars::findLUFourier(ars1.coefficients(), bbbs[i].x0, bbbs[i].x1, bbbs[i].y0, bbbs[i].y1);
        std::cout << i << ": x0 " << RAD2DEG(bbbs[i].x0) << " x1 " << RAD2DEG(bbbs[i].x1) << ", y0 " << bbbs[i].y0 << " y1 " << bbbs[i].y1 << std::endl;
    }


    cuars::FourierOptimizerBB1D optim(ars1.coefficients());
    double xopt, ymin, ymax;
    optim.enableXTolerance(true);
    optim.enableYTolerance(true);
    optim.setXTolerance(M_PI / 180.0 * 0.5);
    optim.setYTolerance(1.0);
    optim.findGlobalMax(0, M_PI, xopt, ymin, ymax);
    std::cout << "\n****\nMaximum in x = " << xopt << " (" << RAD2DEG(xopt) << " deg), maximum between [" << ymin << "," << ymax << "]" << std::endl;

    double xopt2, ymax2;
    cuars::findGlobalMaxBBFourier(ars1.coefficients(), 0, M_PI, M_PI / 180.0 * 0.5, 1.0, xopt2, ymax2);
    std::cout << "  repeated evaluation with findGlobalMaxBBFourier(): maximum in x " << xopt2 << " (" << RAD2DEG(xopt2) << " deg), maximum value " << ymax2 << std::endl;



    //    //Free GPU and CPU memory
    cudaFree(d_coefficientsArs2);
    cudaFree(kernelInput2);
    delete coefficientsArs2;
    //    free(coefficientsArs2); //cpu array
    cudaFree(coeffsMat1);
    cudaFree(kernelInput1);
    cudaFree(d_coeffsArs1);
    delete coeffsArs1;

    return 0;
}

void rangeToPoint(double* ranges, int num, int numPadded, double angleMin, double angleRes, std::vector<cuars::Vec2d>& points) {
    cuars::Vec2d p;
    for (int i = 0; i < numPadded; ++i) {
        if (i < num) {
            double a = angleMin + angleRes * i;
            p.x = ranges[i] * cos(a);
            p.y = ranges[i] * sin(a);
            points.push_back(p);
        } else {
            //padding with zeros

            p.x = 0.0;
            p.y = 0.0;
            points.push_back(p);
        }
        //        std::cout << p.x << " " << p.y << std::endl;
    }
}


