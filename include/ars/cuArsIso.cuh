#include <ars/definitions.h>
#include <ars/ars2d.h>
#include <ars/BBOptimizer1d.h>

#include <thrust/host_vector.h>
#include <thrust/device_malloc.h>

int ceilPow2(int n) {
    ARS_ASSERT(n > 0);

    int exponent = ceil(log2(n));

    int nPadded = std::pow<int>(2, exponent);
    std::cout << "ceilPow2(" << n << ") = " << nPadded << std::endl;



    return nPadded;
}

int sumNaturalsUpToN(int n) {
    ARS_ASSERT(n > 0);

    int result = 0.5 * n * (n + 1);

    std::cout << "sumNaturals(" << n << ") = " << result << std::endl;

    return result;
}

__device__
int getIfromTid(int tid, int n) {
    if (tid < 0 || n < 0)
        return -1;
    const int nId = n - 1; //max n in ids (indices start from 0)
    const int tidStart = tid; //useful for debugging
    /*i is equal to the number of times that we can subtract NID, NID-1, NID-2, ...
     from tid before tid goes below 0*/
    int i = 0;
    if (tid < nId) {
        //        printf("tid %d i %d       n %d\n", tidStart, i, n);
        return i;
    }

    if (tid >= 0.5 * nId * (nId + 1)) { //maybe change this to improve efficiency? 
        //        printf("tid %d i %d       n %d\n", tidStart, i, n);
        i = n;
        return n;
    }

    while (tid >= 0) {
        tid -= (nId - i);
        i++;
    }
    //    printf("tid %d i %d       n %d\n", tidStart, i - 1, n);
    return i - 1;
}

__device__
int getJfromTid(int tid, int n, int i) {
    if (tid < 0 || n < 0 || i < 0)
        return -1;

    const int nId = n - 1; //max n in ids (indices start from 0)
    const int tidStart = tid; //useful for debugging
    const int iStart = i; //useful for debugging
    /*i is equal to the number of times that we can subtract NID, NID-1, NID-2, ...
     from tid before tid goes below 0*/
    int j = -1;

    if (tid >= 0.5 * nId * (nId + 1) || i > nId) { //maybe change this to improve efficiency?
        j = n;
        //        printf("tid %d i %d j %d       n %d\n", tidStart, iStart, j, n);
        return n;
    }

    while (i > 0) {
        tid -= (nId - i);
        i--;
    }
    j = tid + 1;
    //    printf("tid %d i %d j %d       n %d\n", tidStart, iStart, j, n);
    return j;
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

__global__
void iigKernelDownward_old(cuars::Vec2d* means, double sigma1, double sigma2, int numPts, int numPtsAfterPadding, int fourierOrder, int numColsPadded, cuars::ArsKernelIsotropic2d::ComputeMode pnebiMode, double* coeffsMat) {
    //    //    a.insertIsotropicGaussians(points, sigma);
    //
    //    int index = blockIdx.x * blockDim.x + threadIdx.x; //index runs through a single block
    //    int stride = blockDim.x * gridDim.x; //total number of threads in the grid
    //
    //    const int totalNumComparisons = numPtsAfterPadding * numPtsAfterPadding;
    //
    //    for (int tid = index; tid < totalNumComparisons; tid += stride) {
    //
    //        int j = tid % numPtsAfterPadding;
    //        int i = (tid - j) / numPtsAfterPadding;
    //        //        printf("i %d j %d\n", i, j);
    //        //        printf("tid %d i %d j %d tidIJ %d --- numPts %d numPtsAfterPadding %d numColsPadded %d totNumComp %d index %d\n", tid, i, j, i * numPtsAfterPadding + j, numPts, numPtsAfterPadding, numColsPadded, totalNumComparisons, index);
    //
    //        if (i >= numPts || j >= numPts || j <= i)
    //            continue;
    //
    //        cuars::Vec2d vecI = means[i];
    //        cuars::Vec2d vecJ = means[j];
    //
    //        //            isotropicKer_.init(means[i], means[j], sigma);
    //        double dx, dy;
    //        dx = vecJ.x - vecI.x;
    //        dy = vecJ.y - vecI.y;
    //        double phi;
    //
    //        //        if (dx == 0 && dy == 0) {
    //        //                        phi = 0.0; //mathematically undefined
    //        //            //            for (int k = 0; k <= numColsPadded; ++k) {
    //        //            //                int rowIndex = (i * numPtsAfterPadding) + j; //it's more a block index rather than row 
    //        //            //                coeffsMat[rowIndex * numColsPadded + k] = 0.0;
    //        //            //            }
    //        ////            continue;
    //        //
    //        //        } else
    //        phi = atan2(dy, dx);
    //
    //        double sigmaValSq = sigma1 * sigma1 + sigma2 * sigma2;
    //        double lambdaSqNorm = 0.25 * (dx * dx + dy * dy) / sigmaValSq;
    //
    //
    //        //            isotropicKer_.updateFourier(arsfOrder_, coeffs_, w);
    //        double wNorm = 1.0 / (numPts * numPts);
    //        double weight = wNorm / sqrt(2.0 * M_PI * sigmaValSq);
    //
    //
    //
    //        //updating Fourier coefficients (2 modes)
    //        if (pnebiMode == cuars::ArsKernelIsotropic2d::ComputeMode::PNEBI_DOWNWARD) {
    //            //                updateARSF2CoeffRecursDown(lambdaSqNorm, phi, w2, nFourier, coeffs);
    //
    //            double cth2, sth2;
    //            cth2 = cos(2.0 * phi);
    //            sth2 = sin(2.0 * phi);
    //            //                updateARSF2CoeffRecursDown(lambda, cth2, sth2, factor, n, coeffs);
    //
    //
    //
    //
    //            int pnebisSz = fourierOrder + 1;
    //            //TODO: find a better solution instead of hard-coding 21
    //            double pnebis[21]; //Fourier Order + 1
    //            if (pnebis == nullptr)
    //                printf("ERROR ALLOCATING WITH NEW[]!\n");
    //            for (int pn = 0; pn < pnebisSz; ++pn)
    //                pnebis[pn] = 0.0;
    //
    //            double sgn, cth, sth, ctmp, stmp;
    //
    //            // Fourier Coefficients 
    //            //                if (coeffs.size() != 2 * n + 2) {
    //            //                    std::cerr << __FILE__ << "," << __LINE__ << ": invalid size of Fourier coefficients vector " << coeffs.size() << " should be " << (2 * n + 2) << std::endl;
    //            //                    return;
    //            //                }
    //
    //            evaluatePnebiVectorGPU(fourierOrder, lambdaSqNorm, pnebis, pnebisSz);
    //            //                ARS_PRINT(pnebis[0]);
    //
    //            //!!!! factor = w2
    //            double factor = weight;
    //            int rowIndex = (i * numPtsAfterPadding) + j; // = tid
    //            coeffsMat[rowIndex * numColsPadded + 0] += 0.5 * factor * pnebis[0];
    //            //            printf("coeff 0 %f\n", 0.5 * factor * pnebis[0]);
    //
    //
    //            sgn = -1.0;
    //            cth = cth2;
    //            sth = sth2;
    //            //!!!! n in the for below is fourierOrder
    //            for (int k = 1; k <= fourierOrder; ++k) {
    //                //                printf("coeff %d %f\n", 2 * k, factor * pnebis[k] * sgn * cth);
    //                //                printf("coeff %d %f\n", 2 * k + 1, factor * pnebis[k] * sgn * sth);
    //                coeffsMat[(rowIndex * numColsPadded) + (2 * k)] += factor * pnebis[k] * sgn * cth;
    //                coeffsMat[(rowIndex * numColsPadded) + ((2 * k) + 1)] += factor * pnebis[k] * sgn * sth;
    //                sgn = -sgn;
    //                ctmp = cth2 * cth - sth2 * sth;
    //                stmp = sth2 * cth + cth2 * sth;
    //                cth = ctmp;
    //                sth = stmp;
    //            }
    //
    //            delete pnebis;
    //        } else
    //            printf("ERROR: pnebi mode is NOT Downward!\n");
    //
    //
    //
    //    }
}

__global__
void iigKernelDownward(cuars::Vec2d* means, double sigma1, double sigma2, int numPts, int fourierOrder, int numColsPadded, cuars::ArsKernelIsotropic2d::ComputeMode pnebiMode, double* coeffsMat) {
    //    a.insertIsotropicGaussians(points, sigma);

    int index = blockIdx.x * blockDim.x + threadIdx.x; //index runs through a single block
    int stride = blockDim.x * gridDim.x; //total number of threads in the grid

    const int totalNumComparisons = gridDim.x * blockDim.x;

    for (int tid = index; tid < totalNumComparisons; tid += stride) {

        int i = getIfromTid(tid, numPts);
        int j = getJfromTid(tid, numPts, i);
        //        printf("i %d j %d\n", i, j);

        if (i >= numPts || j >= numPts || j <= i) {
            //            printf("BAD INDEXING!!!!\n"); //could actually be ok because of padding... maybe optimization can be further improved in this regard
            continue;
        }

        cuars::Vec2d vecI = means[i];
        cuars::Vec2d vecJ = means[j];

        //            isotropicKer_.init(means[i], means[j], sigma);
        double dx, dy;
        dx = vecJ.x - vecI.x;
        dy = vecJ.y - vecI.y;
        double phi;

        //        if (dx == 0 && dy == 0) {
        //                        phi = 0.0; //mathematically undefined
        //            //            for (int k = 0; k <= numColsPadded; ++k) {
        //            //                int rowIndex = (i * numPtsAfterPadding) + j; //it's more a block index rather than row 
        //            //                coeffsMat[rowIndex * numColsPadded + k] = 0.0;
        //            //            }
        ////            continue;
        //
        //        } else
        phi = atan2(dy, dx);

        double sigmaValSq = sigma1 * sigma1 + sigma2 * sigma2;
        double lambdaSqNorm = 0.25 * (dx * dx + dy * dy) / sigmaValSq;


        //            isotropicKer_.updateFourier(arsfOrder_, coeffs_, w);
        double wNorm = 1.0 / (numPts * numPts);
        double weight = wNorm / sqrt(2.0 * M_PI * sigmaValSq);



        //updating Fourier coefficients (2 modes)
        if (pnebiMode == cuars::ArsKernelIsotropic2d::ComputeMode::PNEBI_DOWNWARD) {
            //                updateARSF2CoeffRecursDown(lambdaSqNorm, phi, w2, nFourier, coeffs);

            double cth2, sth2;
            cth2 = cos(2.0 * phi);
            sth2 = sin(2.0 * phi);
            //                updateARSF2CoeffRecursDown(lambda, cth2, sth2, factor, n, coeffs);




            int pnebisSz = fourierOrder + 1;
            //TODO: find a better solution instead of hard-coding 21
            double pnebis[21]; //Fourier Order + 1
            if (pnebis == nullptr)
                printf("ERROR ALLOCATING WITH NEW[]!\n");
            for (int pn = 0; pn < pnebisSz; ++pn)
                pnebis[pn] = 0.0;

            double sgn, cth, sth, ctmp, stmp;

            // Fourier Coefficients 
            //                if (coeffs.size() != 2 * n + 2) {
            //                    std::cerr << __FILE__ << "," << __LINE__ << ": invalid size of Fourier coefficients vector " << coeffs.size() << " should be " << (2 * n + 2) << std::endl;
            //                    return;
            //                }

            evaluatePnebiVectorGPU(fourierOrder, lambdaSqNorm, pnebis, pnebisSz);
            //                ARS_PRINT(pnebis[0]);

            //!!!! factor = w2
            double factor = weight;
            int rowIndex = tid; // = tid
            coeffsMat[rowIndex * numColsPadded + 0] += 0.5 * factor * pnebis[0];
            //            printf("coeff 0 %f\n", 0.5 * factor * pnebis[0]);


            sgn = -1.0;
            cth = cth2;
            sth = sth2;
            //!!!! n in the for below is fourierOrder
            for (int k = 1; k <= fourierOrder; ++k) {
                //                printf("coeff %d %f\n", 2 * k, factor * pnebis[k] * sgn * cth);
                //                printf("coeff %d %f\n", 2 * k + 1, factor * pnebis[k] * sgn * sth);
                coeffsMat[(rowIndex * numColsPadded) + (2 * k)] += factor * pnebis[k] * sgn * cth;
                coeffsMat[(rowIndex * numColsPadded) + ((2 * k) + 1)] += factor * pnebis[k] * sgn * sth;
                sgn = -sgn;
                ctmp = cth2 * cth - sth2 * sth;
                stmp = sth2 * cth + cth2 * sth;
                cth = ctmp;
                sth = stmp;
            }

            delete pnebis;
        } else
            printf("ERROR: pnebi mode is NOT Downward!\n");



    }
}

__global__
void iigKernelLut(cuars::Vec2d* means, double sigma1, double sigma2, int numPts, int numPtsAfterPadding, int fourierOrder, int numColsPadded, cuars::ArsKernelIsotropic2d::ComputeMode pnebiMode, cuars::PnebiLUT& pnebiLUT, double* coeffsMat) {
    //    a.insertIsotropicGaussians(points, sigma);

    int index = blockIdx.x * blockDim.x + threadIdx.x; //index runs through a single block
    int stride = blockDim.x * gridDim.x; //total number of threads in the grid

    const int totalNumComparisons = numPtsAfterPadding * numPtsAfterPadding;

    for (int tid = index; tid < totalNumComparisons; tid += stride) {

        int j = tid % numPtsAfterPadding;
        int i = (tid - j) / numPtsAfterPadding;
        //        printf("i %d j %d\n", i, j);
        //        printf("tid %d i %d j %d tidIJ %d --- numPts %d numPtsAfterPadding %d numColsPadded %d totNumComp %d index %d\n", tid, i, j, i * numPtsAfterPadding + j, numPts, numPtsAfterPadding, numColsPadded, totalNumComparisons, index);

        if (i >= numPts || j >= numPts || j <= i)
            continue;

        cuars::Vec2d vecI = means[i];
        cuars::Vec2d vecJ = means[j];

        //            isotropicKer_.init(means[i], means[j], sigma);
        double dx, dy;
        dx = vecJ.x - vecI.x;
        dy = vecJ.y - vecI.y;
        double phi;

        //        if (dx == 0 && dy == 0) {
        //                        phi = 0.0; //mathematically undefined
        //            //            for (int k = 0; k <= numColsPadded; ++k) {
        //            //                int rowIndex = (i * numPtsAfterPadding) + j; //it's more a block index rather than row 
        //            //                coeffsMat[rowIndex * numColsPadded + k] = 0.0;
        //            //            }
        ////            continue;
        //
        //        } else
        phi = atan2(dy, dx);

        double sigmaValSq = sigma1 * sigma1 + sigma2 * sigma2;
        double lambdaSqNorm = 0.25 * (dx * dx + dy * dy) / sigmaValSq;
        printf("lambdaSqNorm %f\n", lambdaSqNorm); //just to avoid seeing warning of unused variable when compiling


        //            isotropicKer_.updateFourier(arsfOrder_, coeffs_, w);
        double wNorm = 1.0 / (numPts * numPts);
        double weight = wNorm / sqrt(2.0 * M_PI * sigmaValSq);



        //updating Fourier coefficients (2 modes)
        if (pnebiMode == cuars::ArsKernelIsotropic2d::ComputeMode::PNEBI_LUT) {
            printf("Method not fully implemented!\n");
            continue;

            //                updateARSF2CoeffRecursDownLUT(lambdaSqNorm_, phi_, w2, nFourier, pnebiLut_, coeffs);
            double cth2, sth2;
            //fastCosSin(2.0 * phi, cth2, sth2); //già commentata nell'originale
            cth2 = cos(2.0 * phi);
            sth2 = sin(2.0 * phi);


            int pnebisSz = fourierOrder + 1;
            double *pnebis = new double[pnebisSz];
            double sgn, cth, sth, ctmp, stmp;


            coeffsMat[0] = 0.5 * weight * pnebis[0]; //factor = w2

            sgn = -1.0;
            cth = cth2;
            sth = sth2;
            for (int k = 1; k <= fourierOrder; ++k) {

                coeffsMat[2 * k] = pnebis[k] * weight * sgn * cth;
                coeffsMat[2 * k + 1] = pnebis[k] * weight * sgn * sth;
                sgn = -sgn;
                ctmp = cth2 * cth - sth2 * sth;
                stmp = sth2 * cth + cth2 * sth;
                cth = ctmp;
                sth = stmp;
            }

            delete pnebis;
        } else
            printf("ERROR: pnebi mode is not LUT!\n");

    }
}

__global__
void sumColumns(double* mat, int nrows, int ncols, double* sums) {
    int index = blockIdx.x * blockDim.x + threadIdx.x; //index runs through a single block
    int stride = blockDim.x * gridDim.x; //total number of threads in the grids

    int totalSz = nrows * nrows*ncols; //!! matrix is considered of size (nrows*nrows)*ncols


    for (int idx = index; idx < totalSz; idx += stride) {
        //        int totalIndex = (((i * nrows) + j) * ncols) + k;
        int k = idx % ncols;
        //        int rowIdx = (idx - k) / ncols;
        //        int j = rowIdx % nrows;
        //        int i = (rowIdx - j) / nrows;
        //        printf("i %d j %d k %d rowIdx %d; accessing mat[%d]\n", i, j, k, rowIdx, idx);
        sums[k] += mat[idx];
    }
}

__global__
void sumColumnsNoPadding(double* mat, int nrows, int ncols, double* sums) {
    int index = blockIdx.x * blockDim.x + threadIdx.x; //index runs through a single block
    int stride = blockDim.x * gridDim.x; //total number of threads in the grids

    int totalSz = nrows*ncols; //matrix is considered of size nrows*ncols, with nrows = sumNaturalsUpToN(numPts)


    for (int idx = index; idx < totalSz; idx += stride) {
        //        int totalIndex = (((i * nrows) + j) * ncols) + k;
        int k = idx % ncols;
        //        int rowIdx = (idx - k) / ncols;
        //        int j = rowIdx % nrows;
        //        int i = (rowIdx - j) / nrows;
        //        printf("i %d j %d k %d rowIdx %d; accessing mat[%d]\n", i, j, k, rowIdx, idx);
        sums[k] += mat[idx];
    }
}