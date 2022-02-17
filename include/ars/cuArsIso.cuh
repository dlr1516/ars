


int ceilPow2(int n) {
    ARS_ASSERT(n > 0);

    int exponent = ceil(log2(n));

    int nPadded = std::pow<int>(2, exponent);
    std::cout << "ceilPow2(" << n << ") = " << nPadded << std::endl;



    return nPadded;
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
void iigKernel(cuars::Vec2d* means, double sigma1, double sigma2, int numPts, int numPtsAfterPadding, int fourierOrder, int numColsPadded, cuars::ArsKernelIsotropic2d::ComputeMode pnebiMode, cuars::PnebiLUT& pnebiLUT, double* coeffsMat) {
    //    a.insertIsotropicGaussians(points, sigma);

    int index = blockIdx.x * blockDim.x + threadIdx.x; //index runs through a single block
    int stride = blockDim.x * gridDim.x; //total number of threads in the grid

    const int totalNumComparisons = numPtsAfterPadding * numPtsAfterPadding;

    for (int tid = index; tid < totalNumComparisons; tid += stride) {

        int j = tid % numPtsAfterPadding;
        int i = (tid-j) / numPtsAfterPadding;
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




            int pnebisSz = fourierOrder+1;
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
            int rowIndex = (i * numPtsAfterPadding) + j; // = tid
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
        } else if (pnebiMode == cuars::ArsKernelIsotropic2d::ComputeMode::PNEBI_LUT) {
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
        }



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