
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

#include "ars/Profiler.h"

#include "ars/mpeg7RW.h"
#include "ars/ars2d.cuh"




#define PRINT_DIM(X) std::cout << #X << " rows " << X.rows() << " cols " << X.cols() << std::endl;
#define RAD2DEG(X) (180.0/M_PI*(X))

double mod180(double angle) {
    return (angle - floor(angle / M_PI) * M_PI);
}

struct BoundInterval {
    double x0;
    double x1;
    double y0;
    double y1;
};

int main(int argc, char **argv) {
    cuars::AngularRadonSpectrum2d arsSrc;
    cuars::AngularRadonSpectrum2d arsDst;
    ArsImgTests::PointReaderWriter pointsSrc;
    ArsImgTests::PointReaderWriter pointsDst;


    rofl::ParamMap params;
    std::string filenameCfg;
    std::string filenameSrc;
    std::string filenameDst;
    std::string filenameRot;
    std::string filenameArsSrc;
    std::string filenameArsDst;
    std::string filenameArsRot;
    std::string filenameArsCor;
    std::string filenameCovSrc;
    std::string filenameCovDst;
    int arsOrder;
    double arsSigma, arsThetaToll;
    double rotTrue, rotArs;
    //The variables below are for I/O related functionalities (plot, etc.) that are highly Eigen-based and are present in the CPU-only ArsImgTests...
    //Maybe implement them later
    //    double sampleRes, sampleAng; 
    //    int sampleNum;
    //    bool saveOn;
    //    bool saveCov;


    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, "");
    std::cout << "config filename: " << filenameCfg << std::endl;
    if (filenameCfg != "") {
        params.read(filenameCfg);
    }

    params.read(argc, argv);
    params.getParam<std::string>("src", filenameSrc, "/home/rimlab/Downloads/mpeg7_point_tests/noise000_occl00_rand000/apple-1_xp0686_yp0967_t059_sigma0001_occl000.txt");
    params.getParam<std::string>("dst", filenameDst, "/home/rimlab/Downloads/mpeg7_point_tests/noise000_occl00_rand000/apple-1_xp0749_yn0521_t090_sigma0001_occl000.txt");
    params.getParam<int>("arsOrder", arsOrder, 20);
    params.getParam<double>("arsSigma", arsSigma, 1.0);
    params.getParam<double>("arsTollDeg", arsThetaToll, 1.0);
    arsThetaToll *= M_PI / 180.0;
    //    params.getParam<double>("sampleResDeg", sampleRes, 0.5);
    //    sampleRes *= M_PI / 180.0;
    //    params.getParam<bool>("saveOn", saveOn, false);
    //    params.getParam<bool>("saveCov", saveCov, false);

    std::cout << "\nParameter values:\n";
    params.write(std::cout);
    std::cout << std::endl;




    // Loads files and computes the rotation
    std::cout << "\n*****\nLoading file \"" << filenameSrc << "\"" << std::endl;
    pointsSrc.load(filenameSrc);
    std::cout << "\n*****\nLoading file \"" << filenameDst << "\"" << std::endl;
    pointsDst.load(filenameDst);
    std::cout << "  points src " << pointsSrc.points().size() << ", points dst " << pointsDst.points().size() << std::endl;




    //ARS parameters setting
    arsSrc.setARSFOrder(arsOrder);
    arsDst.setARSFOrder(arsOrder);
    cuars::ArsKernelIsotropic2d::ComputeMode pnebiMode = cuars::ArsKernelIsotropic2d::ComputeMode::PNEBI_DOWNWARD;
    arsSrc.setComputeMode(pnebiMode);
    arsDst.setComputeMode(pnebiMode);


    //Fourier coefficients mega-matrix computation -> parallelization parameters
    int numPts = std::min<int>(pointsSrc.points().size(), pointsDst.points().size()); //the two should normally be equals
    int numPtsAfterPadding = numPts;
    const int gridTotalSize = sumNaturalsUpToN(numPts - 1); //total number of threads in grid Fourier coefficients grid - BEFORE PADDING
    const int blockSize = 256;
    const int numBlocks = floor(gridTotalSize / blockSize) + 1; //number of blocks in grid (each block contains blockSize threads)
    const int gridTotalSizeAfterPadding = blockSize * numBlocks;

    const int coeffsMatNumCols = 2 * arsOrder + 2;
    //    const int coeffsMatNumColsPadded = ceilPow2(coeffsMatNumCols);
    const int coeffsMatNumColsPadded = coeffsMatNumCols;
    const int coeffsMatTotalSz = gridTotalSizeAfterPadding * coeffsMatNumColsPadded; //sumNaturalsUpToN(numPts - 1) * coeffsMatNumColsPadded
    std::cout << "sum parallelization params: " << std::endl
            << " coeffMatNumCols " << coeffsMatNumCols << " coeffsMatTotalSz " << coeffsMatTotalSz << std::endl;

    //Fourier matrix sum -> parallelization parameters
    const int sumBlockSz = 2 * arsOrder + 2;
    const int sumGridSz = blockSize; //unused for now
    std::cout << "Parallelization params:" << std::endl;
    std::cout << "numPts " << numPts << " blockSize " << blockSize << " numBlocks " << numBlocks
            << " gridTotalSize " << gridTotalSize << " gridTotalSizeAP " << gridTotalSizeAfterPadding << std::endl;
    std::cout << "sumSrcBlockSz " << sumBlockSz << " sumGridSz " << sumGridSz << std::endl;


    std::cout << "\n------\n" << std::endl;

    std::cout << "\n\nCalling kernel functions on GPU\n" << std::endl;


    //ARS SRC -> preparation for kernel calls and kernel calls
    cudaEvent_t startSrc, stopSrc; //timing using CUDA events
    cudaEventCreate(&startSrc);
    cudaEventCreate(&stopSrc);
    cuars::Vec2d * kernelInputSrc;
    cudaMalloc((void**) &kernelInputSrc, numPtsAfterPadding * sizeof (cuars::Vec2d));
    cudaMemcpy(kernelInputSrc, pointsSrc.points().data(), numPtsAfterPadding * sizeof (cuars::Vec2d), cudaMemcpyHostToDevice);

    double *d_coeffsMatSrc;
    cudaMalloc((void**) &d_coeffsMatSrc, coeffsMatTotalSz * sizeof (double));
    cudaMemset(d_coeffsMatSrc, 0.0, coeffsMatTotalSz * sizeof (double));
    //    for (int i = 0; i < coeffsMatTotalSz; ++i) {
    //        coeffsMaSrc1[i] = 0.0;
    //    }

    double* d_partsumsSrc;
    cudaMalloc((void**) &d_partsumsSrc, blockSize * coeffsMatNumColsPadded * sizeof (double));
    cudaMemset(d_partsumsSrc, 0.0, blockSize * coeffsMatNumColsPadded * sizeof (double));

    double* d_coeffsArsSrc;
    cudaMalloc((void**) &d_coeffsArsSrc, coeffsMatNumColsPadded * sizeof (double));
    cudaMemset(d_coeffsArsSrc, 0.0, coeffsMatNumColsPadded * sizeof (double));

    cudaEventRecord(startSrc);
    iigKernelDownward << <numBlocks, blockSize >> >(kernelInputSrc, arsSigma, arsSigma, numPts, arsOrder, coeffsMatNumColsPadded, pnebiMode, d_coeffsMatSrc);
    //    sumColumnsNoPadding << <1, sumBlockSz>> >(coeffsMatSrc, gridTotalSizeAfterPadding, coeffsMatNumColsPadded, d_coeffsArsSrc);
    makePartialSums << < coeffsMatNumColsPadded, blockSize >>> (d_coeffsMatSrc, gridTotalSizeAfterPadding, coeffsMatNumColsPadded, d_partsumsSrc);
    sumColumnsPartialSums << <coeffsMatNumColsPadded, 1 >> >(d_partsumsSrc, blockSize, coeffsMatNumColsPadded, d_coeffsArsSrc);
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

    cudaFree(d_coeffsArsSrc);
    cudaFree(d_partsumsSrc);
    cudaFree(d_coeffsMatSrc);
    cudaFree(kernelInputSrc);
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

    double *d_coeffsMatDst; //magari evitare di fare il delete e poi riallocarla è più efficiente (anche se comunque ci sarebbe poi da settare tutto a 0)
    cudaMalloc((void**) &d_coeffsMatDst, coeffsMatTotalSz * sizeof (double));
    cudaMemset(d_coeffsMatDst, 0.0, coeffsMatTotalSz * sizeof (double));
    //    for (int i = 0; i < coeffsMatTotalSz; ++i) {
    //        coeffsMatDst[i] = 0.0;
    //    }

    double* d_partsumsDst;
    cudaMalloc((void**) &d_partsumsDst, blockSize * coeffsMatNumColsPadded * sizeof (double));
    cudaMemset(d_partsumsDst, 0.0, blockSize * coeffsMatNumColsPadded * sizeof (double));

    double* d_coeffsArsDst;
    cudaMalloc((void**) &d_coeffsArsDst, coeffsMatNumColsPadded * sizeof (double));
    cudaMemset(d_coeffsArsDst, 0.0, coeffsMatNumColsPadded * sizeof (double));

    cudaEventRecord(startDst);
    iigKernelDownward << <numBlocks, blockSize >> >(kernelInputDst, arsSigma, arsSigma, numPts, arsOrder, coeffsMatNumColsPadded, pnebiMode, d_coeffsMatDst);
    //    sumColumnsNoPadding << <1, sumBlockSz>> >(coeffsMatDst, gridTotalSizeAfterPadding, coeffsMatNumColsPadded, d_coeffsArsDst);
    makePartialSums << < coeffsMatNumColsPadded, blockSize >>> (d_coeffsMatDst, gridTotalSizeAfterPadding, coeffsMatNumColsPadded, d_partsumsDst);
    sumColumnsPartialSums << <coeffsMatNumColsPadded, 1 >> >(d_partsumsDst, blockSize, coeffsMatNumColsPadded, d_coeffsArsDst);
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

    cudaFree(d_coeffsArsDst);
    cudaFree(d_partsumsDst);
    cudaFree(d_coeffsMatDst);
    cudaFree(kernelInputDst);
    cudaEventDestroy(startDst);
    cudaEventDestroy(stopDst);
    //END OF ARS DST





    //Computation final computations (correlation, ...) on CPU
    std::cout << "\nARS Coefficients:\n";
    std::cout << "Coefficients: Src, Dst, Cor" << std::endl;

    double thetaMax, corrMax, fourierTol;
    fourierTol = 1.0; // TODO: check for a proper tolerance

    std::vector<double> coeffsCor;
    {
        cuars::ScopedTimer("ars.correlation()");
        std::vector<double> tmpSrc;
        tmpSrc.assign(coeffsArsSrc, coeffsArsSrc + coeffsMatNumColsPadded);
        std::vector<double> tmpDst;
        tmpDst.assign(coeffsArsDst, coeffsArsDst + coeffsMatNumColsPadded);
        cuars::computeFourierCorr(tmpSrc, tmpDst, coeffsCor);
        cuars::findGlobalMaxBBFourier(coeffsCor, 0.0, M_PI, arsThetaToll, fourierTol, thetaMax, corrMax);
        rotArs = thetaMax;
    }



    arsSrc.setCoefficients(coeffsArsSrc, coeffsMatNumCols);
    //    for (int i = 0; i < coeffsVectorMaxSz; i++) {
    //        std::cout << "arsSrc - coeff_d[" << i << "] " << d_coeffsMat1[i] << std::endl;
    //    }
    arsDst.setCoefficients(coeffsArsDst, coeffsMatNumCols);
    for (int i = 0; i < arsSrc.coefficients().size() && i < arsDst.coefficients().size(); ++i) {
        std::cout << "\t" << i << " \t" << arsSrc.coefficients().at(i) << " \t" << arsDst.coefficients().at(i) << " \t" << coeffsCor[i] << std::endl;
    }
    std::cout << std::endl;



    // Computes the rotated points,centroid, affine transf matrix between src and dst
    ArsImgTests::PointReaderWriter pointsRot(pointsSrc.points());
    cuars::Vec2d centroidSrc = pointsSrc.computeCentroid();
    cuars::Vec2d centroidDst = pointsDst.computeCentroid();
    cuars::Affine2d rotSrcDst = ArsImgTests::PointReaderWriter::coordToTransform(0.0, 0.0, rotArs);
    //    cuars::Vec2d translSrcDst = centroidDst - rotSrcDst * centroidSrc;
    cuars::Vec2d translSrcDst;
    cuars::vec2diff(translSrcDst, centroidDst, cuars::aff2TimesVec2WRV(rotSrcDst, centroidSrc));
    //    std::cout << "centroidSrc " << centroidSrc.transpose() << "\n"
    //            << "rotSrcDst\n" << rotSrcDst.matrix() << "\n"
    //            << "translation: [" << translSrcDst.transpose() << "] rotation[deg] " << (180.0 / M_PI * rotArs) << "\n";
    std::cout << "centroidSrc " << centroidSrc.x << " \t" << centroidSrc.y << "\n"
            << "centroidDst " << centroidDst.x << " \t" << centroidDst.y << "\n"
            << "rotSrcDst\n" << rotSrcDst << "\n"
            << "translation: [" << translSrcDst.x << " \t" << translSrcDst.y << "] rotation[deg] " << (180.0 / M_PI * rotArs) << "\n";
    pointsRot.applyTransform(translSrcDst.x, translSrcDst.y, rotArs);



    rotTrue = pointsDst.getRotTheta() - pointsSrc.getRotTheta();
    std::cout << "\n***\npointsDst.getrotTheta() [deg]" << (180 / M_PI * pointsDst.getRotTheta())
            << ", pointsSrc.getrotTheta() [deg] " << (180.0 / M_PI * pointsSrc.getRotTheta()) << "\n";
    std::cout << "rotTrue[deg] \t" << (180.0 / M_PI * rotTrue) << " \t" << (180.0 / M_PI * mod180(rotTrue)) << std::endl;
    std::cout << "rotArs[deg] \t" << (180.0 / M_PI * rotArs) << " \t" << (180.0 / M_PI * mod180(rotArs)) << std::endl;

    //Free CPU memory
    delete coeffsArsSrc;
    delete coeffsArsDst;


    return 0;
}





