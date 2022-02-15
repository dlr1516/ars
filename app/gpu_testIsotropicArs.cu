
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

#include <ars/definitions.h>
#include <ars/ars2d.h>
#include <ars/BBOptimizer1d.h>

#include <thrust/host_vector.h>
#include <thrust/device_malloc.h>

#include <chrono>

#include <device_launch_parameters.h>



#define PRINT_DIM(X) std::cout << #X << " rows " << X.rows() << " cols " << X.cols() << std::endl;
#define RAD2DEG(X) (180.0/M_PI*(X))

struct BoundInterval {
    double x0;
    double x1;
    double y0;
    double y1;
};

void rangeToPoint(double* ranges, int num, int numPadded, double angleMin, double angleRes, std::vector<cuars::Vec2d>& points);

int ceilPow2(int n);

//__device__

double evaluatePnebi0Polynom(double x) {
    double t, t2, tinv, val;

    if (x < 0.0) x = -x;
    t = x / 3.75;

    if (t < 1.0) {
        t2 = t*t;
        val = 1.0 + t2 * (3.5156229 + t2 * (3.0899424 + t2 * (1.2067492 + t2 * (0.2659732 + t2 * (0.360768e-1 + t2 * 0.45813e-2)))));
        val = 2.0 * exp(-x) * val;
    } else {
        tinv = 1 / t;
        val = (0.39894228 + tinv * (0.1328592e-1 + tinv * (0.225319e-2 + tinv * (-0.157565e-2 + tinv *
                (0.916281e-2 + tinv * (-0.2057706e-1 + tinv * (0.2635537e-1 + tinv * (-0.1647633e-1 + tinv * 0.392377e-2))))))));
        val = 2.0 * val / sqrt(x);
    }

    return val;
}

//__device__

void evaluatePnebiVectorGPU(int n, double x, double* pnebis, int pnebisSz) {
    double factor, seqPrev, seqCurr, seqNext;
    //    if (pnebis.size() < n + 1) { //questa condizione dovrebbe essere già garantita prima della chiamata di evaluatePnebiVectorGPU
    //        pnebis.resize(n + 1); //ovvero: il questo resizing non dovrebbe essere necessario
    //    }

    if (x < 0.0) x = -x;

    // If x~=0, then BesselI(0,x) = 1.0 and BesselI(k,x) = 0.0 for k > 0.
    // Thus, PNEBI(0,x) = 2.0 and PNEBI(k,x) = 0.0 for k > 0.
    //TODO 9): this if should be done in iigKernel before calling present function
    if (x < 1e-6) {
        pnebis[0] = 2.0;
        for (int i = 1; i < pnebisSz; ++i)
            pnebis[i] = 0.0;
        return;
    }

    // Computes bessel function using back recursion
    factor = 2.0 / x;
    seqPrev = 0.0; // bip
    seqCurr = 1.0; // bi
    seqNext = 0.0; // bim
    for (int k = 2 * (n + (int) sqrt(40.0 * n)); k >= 0; --k) {
        seqNext = seqPrev + factor * k * seqCurr;
        seqPrev = seqCurr;
        seqCurr = seqNext;
        if (k <= n) {
            pnebis[k] = seqPrev;
        }
        // To avoid overflow!
        if (seqCurr > cuars::BIG_NUM) {
            seqPrev *= cuars::SMALL_NUM;
            seqCurr *= cuars::SMALL_NUM;
            for (int i = 0; i < pnebisSz; ++i) {
                pnebis[i] *= cuars::SMALL_NUM;
            }
            //std::cerr << __FILE__ << "," << __LINE__ << ": ANTI-OVERFLOW!" << std::endl;
        }
    }

    double scaleFactor = evaluatePnebi0Polynom(x) / pnebis[0];
    for (int i = 0; i < pnebisSz; ++i) {
        pnebis[i] = scaleFactor * pnebis[i];
    }
}

//__global__

void iigKernel(cuars::Vec2d* means, double sigma1, double sigma2, int numPts, int numPtsAfterPadding, int fourierOrder, int numColsPadded, cuars::ArsKernelIsotropic2d::ComputeMode pnebiMode, cuars::PnebiLUT& pnebiLUT, double* coeffsMat) {
    //    a.insertIsotropicGaussians(points, sigma);

    //TODO 1): 4 righe sotto: DA FARE NEL MAIN PRIMA DI CHIAMARE LA FUNZIONE; vengono fatte una tantum prima del for
    //    if (coeffs_.size() != 2 * arsfOrder_ + 2) {
    //        coeffs_.resize(2 * arsfOrder_ + 2);
    //    }
    //    std::fill(coeffs_.begin(), coeffs_.end(), 0.0);

    //    int index = blockIdx.x * blockDim.x + threadIdx.x; //index runs through a single block
    //    int stride = blockDim.x * gridDim.x; //total number of threads in the grid

    const int totalNumComparisons = numPtsAfterPadding * numPtsAfterPadding;

    for (int tid = 0; tid < totalNumComparisons; ++tid) {

        int j = tid % numPtsAfterPadding;
        int i = (tid - j) / numPtsAfterPadding;
        printf("i %d j %d\n", i, j);
        //        printf("tid %d i %d j %d tidIJ %d --- numPts %d numPtsAfterPadding %d numColsPadded %d totNumComp %d\n", tid, i, j, i * numPtsAfterPadding + j, numPts, numPtsAfterPadding, numColsPadded, totalNumComparisons);

        if (i >= numPts || j >= numPts)
            continue;

        cuars::Vec2d vecI = means[i];
        cuars::Vec2d vecJ = means[j];

        //            isotropicKer_.init(means[i], means[j], sigma);
        double dx, dy;
        dx = vecJ.x - vecI.x;
        dy = vecJ.y - vecI.y;
        double phi;

        if (dx == 0 && dy == 0) {
            //            phi = 0.0; //mathematically undefined
            //            for (int k = 0; k <= numColsPadded; ++k) {
            //                int rowIndex = (i * numPtsAfterPadding) + j; //it's more a block index rather than row 
            //                coeffsMat[rowIndex * numColsPadded + k] = 0.0;
            //            }
            continue;

        } else
            phi = atan2(dy, dx);

        double sigmaValSq = sigma1 * sigma1 + sigma2 * sigma2;
        double lambdaSqNorm = 0.25 * (dx * dx + dy * dy) / sigmaValSq;


        //            isotropicKer_.updateFourier(arsfOrder_, coeffs_, w);
        double weight = 1.0 / (numPts * numPts);
        double w2 = weight / sqrt(2.0 * M_PI * sigmaValSq);

        //TODO 2): TROVARE UNA SOLUZIONE A QUESTO RESIZING (farlo prima di dimensione fissa sufficiente nel main?)
        //            if (coeffs.size() != 2 * nFourier + 2) {
        //                coeffs.resize(2 * nFourier + 2);
        //            }

        //TODO 3): fare questa inizializzazione della LUT nel main
        //            if (pnebiLut_.getOrderMax() < nFourier) {
        //                ARS_ERROR("LUT not initialized to right order. Initialized now.");
        //                pnebiLut_.init(nFourier, 0.0001);
        //            }

        //updating Fourier coefficients (2 modes)
        if (pnebiMode == cuars::ArsKernelIsotropic2d::ComputeMode::PNEBI_DOWNWARD) {
            //                updateARSF2CoeffRecursDown(lambdaSqNorm, phi, w2, nFourier, coeffs);

            double cth2, sth2;
            cth2 = cos(2.0 * phi);
            sth2 = sin(2.0 * phi);
            //                updateARSF2CoeffRecursDown(lambda, cth2, sth2, factor, n, coeffs);


            //TODO 4): make pnebis a double*
            //can solve it with cuda/gpu malloc here? otherwise just pass the needed pointer to the function?
            //for now I just declare it here
            //                std::vector<double> pnebis(n + 1);
            int pnebisSz = fourierOrder + 1;
            double *pnebis = new double[pnebisSz];

            double sgn, cth, sth, ctmp, stmp;

            // Fourier Coefficients 
            //                if (coeffs.size() != 2 * n + 2) {
            //                    std::cerr << __FILE__ << "," << __LINE__ << ": invalid size of Fourier coefficients vector " << coeffs.size() << " should be " << (2 * n + 2) << std::endl;
            //                    return;
            //                }

            //                                TODO 5): expand evaluatePnebiVector() below
            evaluatePnebiVectorGPU(fourierOrder, lambdaSqNorm, pnebis, pnebisSz);
            //                ARS_PRINT(pnebis[0]);

            //!!!! factor = w2
            double factor = w2;
            int rowIndex = (i * numPtsAfterPadding) + j; // = tid
            coeffsMat[rowIndex * numColsPadded + 0] += 0.5 * factor * pnebis[0];
            std::cout << "coeff0" << 0.5 * factor * pnebis[0] << std::endl;


            sgn = -1.0;
            cth = cth2;
            sth = sth2;
            //!!!! n in the for below is fourierOrder
            //                for (int k = 1; k <= n; ++k) {
            for (int k = 1; k <= fourierOrder; ++k) {
                std::cout << "coeff" << 2 * k << " " << factor * pnebis[k] * sgn * cth << std::endl;
                std::cout << "coeff" << 2 * k + 1 << " " << factor * pnebis[k] * sgn * cth << std::endl;
                coeffsMat[(rowIndex * numColsPadded) + (2 * k)] += factor * pnebis[k] * sgn * cth;
                coeffsMat[(rowIndex * numColsPadded) + ((2 * k) + 1)] += factor * pnebis[k] * sgn * sth;
                sgn = -sgn;
                ctmp = cth2 * cth - sth2 * sth;
                stmp = sth2 * cth + cth2 * sth;
                cth = ctmp;
                sth = stmp;
            }
        } else if (pnebiMode == cuars::ArsKernelIsotropic2d::ComputeMode::PNEBI_LUT) {
            //                updateARSF2CoeffRecursDownLUT(lambdaSqNorm_, phi_, w2, nFourier, pnebiLut_, coeffs);
            double cth2, sth2;
            //fastCosSin(2.0 * phi, cth2, sth2); //già commentata nell'originale
            cth2 = cos(2.0 * phi);
            sth2 = sin(2.0 * phi);

            //TODO 6): find a workaround for this pnebis vector; for now I just initialize here a double* pnebis;
            //                std::vector<double> pnebis(fourierOrder + 1); //prima riga della funzione omonima chiamata da dentro l'inline
            int pnebisSz = fourierOrder + 1;
            double *pnebis = new double[pnebisSz];
            double sgn, cth, sth, ctmp, stmp;

            //TODO 7): minor problem... seems just to be a check of standing conditions. Still... might be useful to understand it in order to fix dimensions of pointers passed to iigKernel
            // Fourier Coefficients 
            //                if (coeffs.size() != 2 * fourierOrder + 2 || pnebiLUT.getOrderMax() < fourierOrder) {
            //                    std::cerr << __FILE__ << "," << __LINE__ << ": one of these conditions failed:"
            //                            << "\n  size of Fourier coefficients vector " << coeffs.size() << " should be " << (2 * n + 2)
            //                            << "\n  LUT max order is " << pnebiLUT.getOrderMax() << " >= " << n
            //                            << std::endl;
            //                    return;
            //                }

            // TODO 8): SOLVE PROBLEM OF FUNCTION COMMENTED BELOW (NOTE THAT ITS CODE HAS ALREADY BEEN COPIED IN THE SCOPE BELOW THE COMMENTED CALLING OF THE FUNCTION)
            //                pnebiLUT.eval(lambdaSqNorm, pnebis);
            //                evalPnebiLUT2();
            //                //ARS_PRINT(pnebis[0]);


            coeffsMat[0] = 0.5 * w2 * pnebis[0]; //factor = w2

            sgn = -1.0;
            cth = cth2;
            sth = sth2;
            for (int k = 1; k <= fourierOrder; ++k) {

                coeffsMat[2 * k] = pnebis[k] * w2 * sgn * cth;
                coeffsMat[2 * k + 1] = pnebis[k] * w2 * sgn * sth;
                sgn = -sgn;
                ctmp = cth2 * cth - sth2 * sth;
                stmp = sth2 * cth + cth2 * sth;
                cth = ctmp;
                sth = stmp;
            }
        }



    }
}

int main(void) {
    double acesRanges[] =  {4.32, 1.13, 3.51, 2.54, 4.25, 2.17, 4.85, 1.27, 7.24, 9.43, 1.36, 6.30, 8.36, 7.61, 0.31, 7.49, 0.38, 7.23, 3.97, 0.54, 1.38, 9.44, 5.93, 7.57, 3.96, 8.19, 1.44, 1.73, 8.01, 3.85, 4.58, 3.71, 2.28, 3.79, 5.43, 3.57, 9.24, 8.47, 4.52, 0.60, 1.07, 0.56, 3.26, 8.90, 0.95, 4.48, 8.78, 1.41, 4.63, 0.64, 4.19, 6.31, 9.86, 1.68, 9.03, 6.51, 8.70, 7.58, 0.10, 1.35, 3.06, 1.72, 5.98, 3.66, 1.18, 5.54, 8.98, 3.52, 6.17, 7.10, 6.26, 4.23, 6.18, 2.06, 4.27, 2.21, 7.52, 6.30, 8.71, 3.17, 8.56, 4.65, 0.16, 2.02, 7.05, 6.34, 6.37, 0.66, 4.33, 1.10, 9.50, 4.68, 4.72, 4.55, 0.69, 7.38, 3.77, 3.22, 1.43, 9.43, 6.41, 5.53, 7.00, 4.61, 5.42, 3.80, 1.73, 5.78, 0.45, 6.42, 4.99, 9.01, 8.86, 1.86, 4.33, 4.29, 3.00, 3.60, 4.66, 7.72, 5.54, 0.75, 2.48, 2.14, 7.02, 2.14, 1.10, 4.36, 4.06, 6.81, 1.24, 7.58, 9.29, 9.41, 0.83, 9.31, 0.24, 1.57, 6.17, 9.61, 3.38, 9.74, 6.89, 2.30, 5.51, 9.17, 5.11, 8.25, 0.72, 0.17, 8.61, 7.46, 6.39, 1.24, 7.01, 3.94, 0.08, 4.82, 3.86, 1.05, 3.05, 5.37, 6.21, 2.88, 6.86, 6.00, 1.17, 8.85, 0.23, 2.72, 7.51, 8.84, 6.77, 7.18, 1.79, 8.43, 3.02, 5.86, 5.36, 9.83};
    cuars::AngularRadonSpectrum2d ars1;
    cuars::AngularRadonSpectrum2d ars2;
    std::chrono::system_clock::time_point timeStart, timeStop;
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

    //conversion
    std::vector<cuars::Vec2d> acesPointsSTL;
    rangeToPoint(acesRanges, numPts, numPtsAfterPadding, -0.5 * M_PI, M_PI / 180.0 * 1.0, acesPointsSTL);

    thrust::host_vector<cuars::Vec2d> acesPointsHost(acesPointsSTL.begin(), acesPointsSTL.end());


    //    cuars::Vec2d firstElement; //??



    timeStart = std::chrono::system_clock::now();
    //ars1 kernel call
    //    ars1.insertIsotropicGaussians(acesPoints1, sigma);
    cuars::Vec2d * kernelInput1 = new cuars::Vec2d [numPtsAfterPadding];
    //    cudaMallocManaged((void**) &kernelInput1, numPtsAfterPadding * sizeof (cuars::Vec2d));
    //    cudaMemcpy(kernelInput1, acesPointsHost.data(), numPtsAfterPadding * sizeof (cuars::Vec2d), cudaMemcpyDefault);
    for (int i = 0; i < numPtsAfterPadding; ++i) {
        kernelInput1[i] = acesPointsSTL[i];
    }

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
    double *coeffsMat1 = new double [coeffsMatTotalSz];
    //    cudaMallocManaged((void**) &coeffsMat1, coeffsMatTotalSz * sizeof (double));
    //    cudaMemset(coeffsMat1, 0.0, coeffsMatTotalSz * sizeof (double));
    for (int i = 0; i < coeffsMatTotalSz; ++i) {
        coeffsMat1[i] = 0.0;
    }


    cuars::PnebiLUT pnebiLUT1; //LUT setup
    double lutPrecision = 0.001; //LUT setup
    pnebiLUT1.init(fourierOrder, lutPrecision); //LUT setup
    if (pnebiLUT1.getOrderMax() < fourierOrder) { //LUT setup
        ARS_ERROR("LUT not initialized to right order. Initialized now."); //LUT setup
        pnebiLUT1.init(fourierOrder, 0.0001); //LUT setup
    }


    //    iigKernel << < 1, 1 >> >(kernelInput1, sigma, sigma, numPts, numPtsAfterPadding, fourierOrder, coeffsMatNumColsPadded, pnebiMode, pnebiLUT1, coeffsMat1);
    iigKernel(kernelInput1, sigma, sigma, numPts, numPtsAfterPadding, fourierOrder, coeffsMatNumColsPadded, pnebiMode, pnebiLUT1, coeffsMat1);

    //    cudaMemcpy(coeffsMat1, d_coeffsMat1, coeffsVectorMaxSz * sizeof (double), cudaMemcpyDefault);
    //end of kernel call


    double* coeffsArs1 = new double [coeffsMatNumColsPadded];
    for (int k = 0; k < coeffsMatNumColsPadded; ++k)
        coeffsArs1[k] = 0.0; //init coeffsArs vector to 0    
    for (int i = 0; i < numPtsAfterPadding; ++i)
        for (int j = 0; j < numPtsAfterPadding; ++j)
            for (int k = 0; k < coeffsMatNumColsPadded; ++k) {
                int totalIndex = (((i * numPtsAfterPadding) + j) * coeffsMatNumColsPadded) + k;
                coeffsArs1[k] += coeffsMat1[totalIndex];
            }
    for (int i = 0; i < coeffsMatNumColsPadded; ++i) {
        std::cout << "coeffsArs1[" << i << "] " << coeffsArs1[i] << std::endl;
    }

    timeStop = std::chrono::system_clock::now();
    double timeArs1 = (double) std::chrono::duration_cast<std::chrono::milliseconds>(timeStop - timeStart).count();
    cudaDeviceSynchronize();
    std::cout << "insertIsotropicGaussians() " << timeArs1 << " ms" << std::endl;
    //END OF ARS1

    std::cout << "\n------\n" << std::endl;


    //ARS2    
    ars2.setComputeMode(cuars::ArsKernelIsotropic2d::ComputeMode::PNEBI_LUT);

    timeStart = std::chrono::system_clock::now();

    //kernel call
    //    ars2.insertIsotropicGaussians(acesPoints1, sigma);
    cuars::Vec2d* kernelInput2;
    cudaMalloc((void **) &kernelInput2, numPtsAfterPadding * sizeof (cuars::Vec2d));
    pnebiMode = cuars::ArsKernelIsotropic2d::ComputeMode::PNEBI_LUT;

    double *coefficientsArs2 = new double[coeffsMatTotalSz](); //() initialize to 0
    double *d_coefficientsArs2; //d_ stands for device
    //    const int coeffsVectorMaxSz = 2 * fourierOrder + 2; //already initialized in ars1
    cudaMalloc(&d_coefficientsArs2, coeffsMatTotalSz * sizeof (double)); //maybe directly use cudaMemset?
    cudaMemcpy(d_coefficientsArs2, coefficientsArs2, coeffsMatTotalSz * sizeof (double), cudaMemcpyHostToDevice);

    cuars::PnebiLUT pnebiLUT2; //LUT setup
    //    double lutPrecision = 0.001; //already initialized for pnebiLUT1
    pnebiLUT2.init(fourierOrder, lutPrecision); //LUT setup
    if (pnebiLUT2.getOrderMax() < fourierOrder) { //LUT setup
        ARS_ERROR("LUT not initialized to right order. Initialized now."); //LUT setup
        pnebiLUT2.init(fourierOrder, 0.0001); //LUT setup
    }

    //    iigKernel << < numBlocks, blockSize >> >(thrust::raw_pointer_cast<ars::Vec2d*>(kernelInput2.data()), sigma, sigma, numPts, paddedPtVecSz, fourierOrder, pnebiMode, pnebiLUT2, d_coefficientsArs2);
    cudaMemcpy(coefficientsArs2, d_coefficientsArs2, coeffsMatTotalSz * sizeof (double), cudaMemcpyDeviceToHost);
    //end of kernel call for ARS2



    timeStop = std::chrono::system_clock::now();
    double timeArs2 = (double) std::chrono::duration_cast<std::chrono::milliseconds>(timeStop - timeStart).count();
    cudaDeviceSynchronize();
    std::cout << "insertIsotropicGaussians() " << timeArs2 << " ms" << std::endl;
    //END OF ARS1




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
    //    cudaFree(d_coefficientsArs2);
    //    cudaFree(kernelInput2);
    //    //    free(coefficientsArs2); //array
    //    cudaFree(coeffsMat1);
    //    cudaFree(kernelInput1);
    //    //    free(coeffsArs1);

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

int ceilPow2(int n) {
    ARS_ASSERT(n > 0);

    int exponent = ceil(log2(n));

    int nPadded = std::pow<int>(2, exponent);
    std::cout << "Number of points: " << n << " -> afeer padding = " << nPadded << std::endl;



    return nPadded;
}


