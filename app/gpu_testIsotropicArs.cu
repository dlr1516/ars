
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

#include <thrust/device_vector.h>

#include <chrono>
#include <sys/param.h>
#include <device_launch_parameters.h>



#define PRINT_DIM(X) std::cout << #X << " rows " << X.rows() << " cols " << X.cols() << std::endl;
#define RAD2DEG(X) (180.0/M_PI*(X))

struct BoundInterval {
    double x0;
    double x1;
    double y0;
    double y1;
};

double acesRanges[] = {50.00, 50.00, 50.00, 5.26, 5.21, 5.06, 5.01, 3.01, 2.94, 2.89, 2.84, 2.74, 2.69, 2.64, 2.59, 2.54, 2.49, 2.49, 2.44, 2.39, 2.34, 2.29, 2.29, 2.29, 2.39, 2.39, 2.49, 2.51, 2.61, 2.66, 2.76, 2.81, 2.96, 3.01, 3.11, 3.26, 3.01, 3.01, 3.01, 3.06, 3.21, 6.86, 6.86, 6.81, 6.76, 6.71, 6.71, 6.66, 6.61, 6.66, 6.56, 6.56, 6.56, 6.46, 6.46, 6.41, 6.46, 6.46, 4.11, 3.96, 3.96, 4.96, 4.86, 5.21, 7.41, 4.61, 5.16, 6.26, 6.26, 6.31, 4.86, 5.01, 5.86, 5.81, 4.21, 4.26, 4.31, 4.41, 4.39, 4.46, 5.31, 5.06, 5.26, 4.96, 6.01, 5.76, 5.61, 5.36, 5.26, 5.01, 4.21, 4.16, 4.01, 3.91, 3.61, 3.21, 3.26, 3.16, 3.06, 3.01, 3.31, 3.21, 3.16, 2.16, 2.19, 2.16, 2.21, 2.11, 2.01, 2.01, 2.06, 2.84, 2.91, 2.91, 3.01, 3.11, 3.21, 3.81, 4.06, 7.11, 7.06, 7.01, 6.96, 6.86, 4.31, 6.76, 6.71, 6.66, 6.61, 5.46, 5.41, 6.46, 6.21, 6.31, 6.51, 7.26, 7.46, 50.00, 2.01, 1.94, 1.94, 1.94, 2.31, 1.86, 1.84, 1.84, 1.81, 1.96, 26.46, 20.76, 2.11, 2.12, 2.17, 2.14, 2.09, 2.09, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.19, 2.19, 2.24, 2.24, 2.24, 2.24, 2.29, 2.29, 2.29, 2.29, 2.29, 2.39, 2.39, 2.39, 2.44};


void rangeToPoint(double* ranges, int num, int numPadded, double angleMin, double angleRes, thrust::device_vector<ars::Vec2d>& points);

int ceilPow2(int n) {
    ARS_ASSERT(n > 0);

    int exponent = ceil(sqrt(n));

    return (int) std::pow<int>(2, exponent);
}

__device__
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

__device__
void evaluatePnebiVectorGPU(int n, double x, double* pnebis, int pnebisSz) {
    double factor, seqPrev, seqCurr, seqNext;
    //    if (pnebis.size() < n + 1) { //questa condizione dovrebbe essere già garantita prima della chiamata di evaluatePnebiVectorGPU
    //        pnebis.resize(n + 1); //ovvero: il questo resizing non dovrebbe essere necessario
    //    }

    if (x < 0.0) x = -x;

    // If x~=0, then BesselI(0,x) = 1.0 and BesselI(k,x) = 0.0 for k > 0.
    // Thus, PNEBI(0,x) = 2.0 and PNEBI(k,x) = 0.0 for k > 0.
    //TODO 9): this if should be done in iigKernel before calling present function
    //    if (x < 1e-6) {
    //        std::fill(pnebis.begin(), pnebis.end(), 0.0);
    //        pnebis[0] = 2.0;
    //        return;
    //    }

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
        if (seqCurr > ars::BIG_NUM) {
            seqPrev *= ars::SMALL_NUM;
            seqCurr *= ars::SMALL_NUM;
            for (int i = 0; i < pnebisSz; ++i) {
                pnebis[i] *= ars::SMALL_NUM;
            }
            //std::cerr << __FILE__ << "," << __LINE__ << ": ANTI-OVERFLOW!" << std::endl;
        }
    }

    double scaleFactor = evaluatePnebi0Polynom(x) / pnebis[0];
    for (int i = 0; i < pnebisSz; ++i) {
        pnebis[i] = scaleFactor * pnebis[i];
    }
}

__global__
void iigKernel(ars::Vec2d* means, double sigma1, double sigma2, size_t numPtsAfterPadding, int fourierOrder, ars::ArsKernelIsotropic2d::ComputeMode pnebiMode, ars::PnebiLUT& pnebiLUT, double* coefficients) {
    //    a.insertIsotropicGaussians(points, sigma);

    //TODO 1): 4 righe sotto: DA FARE NEL MAIN PRIMA DI CHIAMARE LA FUNZIONE; vengono fatte una tantum prima del for
    //    if (coeffs_.size() != 2 * arsfOrder_ + 2) {
    //        coeffs_.resize(2 * arsfOrder_ + 2);
    //    }
    //    std::fill(coeffs_.begin(), coeffs_.end(), 0.0);

    size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    size_t stride = 8 * 32;


    //    size_t index = blockIdx.x * blockDim.x + threadIdx.x; //index runs through a single block
    //    size_t stride = blockDim.x * gridDim.x; // = total number of threads in the grid
    if (tid < numPtsAfterPadding * numPtsAfterPadding) {
        //    for (size_t totId = tid; totId < 256*256; totId += stride) {

        size_t j = tid % numPtsAfterPadding;
        size_t i = (tid - j) / numPtsAfterPadding; // equivalent of floor(totId/numPtsPadded)
        ars::Vec2d vecI = means[i];
        ars::Vec2d vecJ = means[j];

        //            isotropicKer_.init(means[i], means[j], sigma);
        double dx, dy;
        dx = vecJ.x - vecI.x;
        dy = vecJ.y - vecI.y;
        double phi = atan2(dy, dx);
        double sigmaValSq = sigma1 * sigma1 + sigma2 * sigma2;
        double lambdaSqNorm = 0.25 * (dx * dx + dy * dy) / sigmaValSq;


        //            isotropicKer_.updateFourier(arsfOrder_, coeffs_, w);
        double weight = 1.0 / (numPtsAfterPadding * numPtsAfterPadding);
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

        if (pnebiMode == ars::ArsKernelIsotropic2d::ComputeMode::PNEBI_DOWNWARD) {
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
            coefficients[0] += 0.5 * factor * pnebis[0];
            sgn = -1.0;
            cth = cth2;
            sth = sth2;
            //!!!! n in the for below is fourierOrder
            //                for (int k = 1; k <= n; ++k) {
            for (int k = 1; k <= fourierOrder; ++k) {

                coefficients[2 * k] += factor * pnebis[k] * sgn * cth;
                coefficients[2 * k + 1] += factor * pnebis[k] * sgn * sth;
                sgn = -sgn;
                ctmp = cth2 * cth - sth2 * sth;
                stmp = sth2 * cth + cth2 * sth;
                cth = ctmp;
                sth = stmp;
            }
        } else if (pnebiMode == ars::ArsKernelIsotropic2d::ComputeMode::PNEBI_LUT) {
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


            coefficients[0] += 0.5 * w2 * pnebis[0]; //factor = w2
            sgn = -1.0;
            cth = cth2;
            sth = sth2;
            for (int k = 1; k <= fourierOrder; ++k) {
                coefficients[2 * k] += pnebis[k] * w2 * sgn * cth;
                coefficients[2 * k + 1] += pnebis[k] * w2 * sgn * sth;
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
    ars::AngularRadonSpectrum2d ars1;
    ars::AngularRadonSpectrum2d ars2;
    thrust::device_vector<ars::Vec2d> acesPoints;
    std::chrono::system_clock::time_point timeStart, timeStop;
    double sigma = 0.05;
    int fourierOrder = 20;

    ars1.setARSFOrder(fourierOrder);
    ars2.setARSFOrder(fourierOrder);

    //parallelization parameters
    size_t numPts = 180; // = acesRanges.size()
    const size_t paddedPtVecSz = ceilPow2(numPts);
    const size_t blockSize = 32;
    const size_t numBlocks = (paddedPtVecSz * paddedPtVecSz) / blockSize;
    const size_t gridTotalSize = blockSize*numBlocks;
    rangeToPoint(acesRanges, numPts, paddedPtVecSz, -0.5 * M_PI, M_PI / 180.0 * 1.0, acesPoints);
    //        acesPoints1.push_back(ars::Vector2::Zero());
    ars::Vec2d firstElement;
    //    firstElement.x = 0.0;
    //    firstElement.y = 0.0;
    //    acesPoints.push_back(firstElement);
    std::cout << "Number of input points: " << acesPoints.size() << std::endl;
    //    for (int i = 0; i < acesPoints1.size(); ++i) {
    //        std::cout << i << "\t" << acesPoints1[i].x() << "\t" << acesPoints1[i].y() << std::endl;
    //    }

    //    const int trialNum = 1;
    ars1.initLUT(0.0001);
    ars1.setComputeMode(ars::ArsKernelIsotropic2d::ComputeMode::PNEBI_DOWNWARD);



    timeStart = std::chrono::system_clock::now();

    //ars1 kernel call
    //    ars1.insertIsotropicGaussians(acesPoints1, sigma);
    ars::ArsKernelIsotropic2d::ComputeMode pnebiMode = ars::ArsKernelIsotropic2d::ComputeMode::PNEBI_DOWNWARD;


    ars::VecVec2d kernelInput1(paddedPtVecSz);


    const size_t coeffsVectorMaxSz = 2 * fourierOrder + 2;
    double *coefficientsArs1 = new double[coeffsVectorMaxSz](); //() initialize to 0
    double *d_coefficientsArs1; //d_ stands for device
    cudaMalloc(&d_coefficientsArs1, coeffsVectorMaxSz * sizeof (double));
    cudaMemcpy(d_coefficientsArs1, coefficientsArs1, coeffsVectorMaxSz * sizeof (double), cudaMemcpyHostToDevice);

    ars::PnebiLUT pnebiLUT1; //LUT setup
    double lutPrecision = 0.001; //LUT setup
    pnebiLUT1.init(fourierOrder, lutPrecision); //LUT setup
    if (pnebiLUT1.getOrderMax() < fourierOrder) { //LUT setup
        ARS_ERROR("LUT not initialized to right order. Initialized now."); //LUT setup
        pnebiLUT1.init(fourierOrder, 0.0001); //LUT setup
    }


    iigKernel << < numBlocks, blockSize >> >(thrust::raw_pointer_cast<ars::Vec2d*>(kernelInput1.data()), sigma, sigma, paddedPtVecSz, fourierOrder, pnebiMode, pnebiLUT1, d_coefficientsArs1);
    cudaMemcpy(coefficientsArs1, d_coefficientsArs1, coeffsVectorMaxSz * sizeof (double), cudaMemcpyDeviceToHost);
    //end of kernel call

    timeStop = std::chrono::system_clock::now();
    double timeArs1 = (double) std::chrono::duration_cast<std::chrono::milliseconds>(timeStop - timeStart).count();
    cudaDeviceSynchronize();
    std::cout << "insertIsotropicGaussians() " << timeArs1 << " ms" << std::endl;


    //END OF ARS1


    std::cout << "\n------\n" << std::endl;


    //ARS2    
    ars2.setComputeMode(ars::ArsKernelIsotropic2d::ComputeMode::PNEBI_LUT);

    timeStart = std::chrono::system_clock::now();

    //kernel call
    //    ars2.insertIsotropicGaussians(acesPoints1, sigma);
    ars::Vec2d* kernelInput2 = thrust::raw_pointer_cast(acesPoints.data());
    numPts = acesPoints.size(); //for this dummy example atleast
    pnebiMode = ars::ArsKernelIsotropic2d::ComputeMode::PNEBI_LUT;

    double *coefficientsArs2 = new double[coeffsVectorMaxSz](); //() initialize to 0
    double *d_coefficientsArs2; //d_ stands for device
    //    const size_t coeffsVectorMaxSz = 2 * fourierOrder + 2; //already initialized in ars1
    cudaMalloc(&d_coefficientsArs2, coeffsVectorMaxSz * sizeof (double)); //maybe directly use cudaMemset?
    cudaMemcpy(d_coefficientsArs2, coefficientsArs2, coeffsVectorMaxSz * sizeof (double), cudaMemcpyHostToDevice);

    ars::PnebiLUT pnebiLUT2; //LUT setup
    //    double lutPrecision = 0.001; //already initialized for pnebiLUT1
    pnebiLUT2.init(fourierOrder, lutPrecision); //LUT setup
    if (pnebiLUT2.getOrderMax() < fourierOrder) { //LUT setup
        ARS_ERROR("LUT not initialized to right order. Initialized now."); //LUT setup
        pnebiLUT2.init(fourierOrder, 0.0001); //LUT setup
    }

    //    iigKernel << < numBlocks, blockSize >> >(thrust::raw_pointer_cast<ars::Vec2d*>(kernelInput2.data()), sigma, sigma, paddedPtVecSz, fourierOrder, pnebiMode, pnebiLUT2, d_coefficientsArs2);
    cudaMemcpy(coefficientsArs2, d_coefficientsArs2, coeffsVectorMaxSz * sizeof (double), cudaMemcpyDeviceToHost);
    //end of kernel call for ARS2


    timeStop = std::chrono::system_clock::now();
    double timeArs2 = (double) std::chrono::duration_cast<std::chrono::milliseconds>(timeStop - timeStart).count();
    std::cout << "insertIsotropicGaussians() " << timeArs2 << " ms" << std::endl;



    std::cout << "\nARS Coefficients:\n";
    std::cout << "\ti \tDownward \tLUT\n";
    ars1.setCoefficients(coefficientsArs1, coeffsVectorMaxSz);
    ars2.setCoefficients(coefficientsArs2, coeffsVectorMaxSz);
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
        ars::findLUFourier(ars1.coefficients(), bbbs[i].x0, bbbs[i].x1, bbbs[i].y0, bbbs[i].y1);
        std::cout << i << ": x0 " << RAD2DEG(bbbs[i].x0) << " x1 " << RAD2DEG(bbbs[i].x1) << ", y0 " << bbbs[i].y0 << " y1 " << bbbs[i].y1 << std::endl;
    }


    ars::FourierOptimizerBB1D optim(ars1.coefficients());
    double xopt, ymin, ymax;
    optim.enableXTolerance(true);
    optim.enableYTolerance(true);
    optim.setXTolerance(M_PI / 180.0 * 0.5);
    optim.setYTolerance(1.0);
    optim.findGlobalMax(0, M_PI, xopt, ymin, ymax);
    std::cout << "\n****\nMaximum in x = " << xopt << " (" << RAD2DEG(xopt) << " deg), maximum between [" << ymin << "," << ymax << "]" << std::endl;

    double xopt2, ymax2;
    ars::findGlobalMaxBBFourier(ars1.coefficients(), 0, M_PI, M_PI / 180.0 * 0.5, 1.0, xopt2, ymax2);
    std::cout << "  repeated evaluation with findGlobalMaxBBFourier(): maximum in x " << xopt2 << " (" << RAD2DEG(xopt2) << " deg), maximum value " << ymax2 << std::endl;



    //Free GPU and CPU memory
    cudaFree(d_coefficientsArs2);
    //    free(coefficientsArs2); //array
    cudaFree(d_coefficientsArs1);
    //    free(coefficientsArs1);

    return 0;
}

void rangeToPoint(double* ranges, int num, int numPadded, double angleMin, double angleRes, thrust::device_vector<ars::Vec2d>& points) {
    ars::Vec2d p;
    for (int i = 0; i < numPadded; ++i) {
        if (i < num) {
            double a = angleMin + angleRes * i;
            p.x = ranges[i] * cos(a);
            p.y = ranges[i] * sin(a);
            points.push_back(p);
            //                std::cout << p.x << " " << p.y << std::endl;
        } else {
            //padding with zeros
            p.x = 0.0;
            p.y = 0.0;
            points.push_back(p);
        }

    }
}



