
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


#define PRINT_DIM(X) std::cout << #X << " rows " << X.rows() << " cols " << X.cols() << std::endl;
#define RAD2DEG(X) (180.0/M_PI*(X))

struct BoundInterval {
    double x0;
    double x1;
    double y0;
    double y1;
};

double acesRanges1[] = {50.00, 50.00, 50.00, 5.26, 5.21, 5.06, 5.01, 3.01, 2.94, 2.89, 2.84, 2.74, 2.69, 2.64, 2.59, 2.54, 2.49, 2.49, 2.44, 2.39, 2.34, 2.29, 2.29, 2.29, 2.39, 2.39, 2.49, 2.51, 2.61, 2.66, 2.76, 2.81, 2.96, 3.01, 3.11, 3.26, 3.01, 3.01, 3.01, 3.06, 3.21, 6.86, 6.86, 6.81, 6.76, 6.71, 6.71, 6.66, 6.61, 6.66, 6.56, 6.56, 6.56, 6.46, 6.46, 6.41, 6.46, 6.46, 4.11, 3.96, 3.96, 4.96, 4.86, 5.21, 7.41, 4.61, 5.16, 6.26, 6.26, 6.31, 4.86, 5.01, 5.86, 5.81, 4.21, 4.26, 4.31, 4.41, 4.39, 4.46, 5.31, 5.06, 5.26, 4.96, 6.01, 5.76, 5.61, 5.36, 5.26, 5.01, 4.21, 4.16, 4.01, 3.91, 3.61, 3.21, 3.26, 3.16, 3.06, 3.01, 3.31, 3.21, 3.16, 2.16, 2.19, 2.16, 2.21, 2.11, 2.01, 2.01, 2.06, 2.84, 2.91, 2.91, 3.01, 3.11, 3.21, 3.81, 4.06, 7.11, 7.06, 7.01, 6.96, 6.86, 4.31, 6.76, 6.71, 6.66, 6.61, 5.46, 5.41, 6.46, 6.21, 6.31, 6.51, 7.26, 7.46, 50.00, 2.01, 1.94, 1.94, 1.94, 2.31, 1.86, 1.84, 1.84, 1.81, 1.96, 26.46, 20.76, 2.11, 2.12, 2.17, 2.14, 2.09, 2.09, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.19, 2.19, 2.24, 2.24, 2.24, 2.24, 2.29, 2.29, 2.29, 2.29, 2.29, 2.39, 2.39, 2.39, 2.44};


void rangeToPoint(double* ranges, int num, double angleMin, double angleRes, thrust::device_vector<ars::Vec2d>& points);

__global__
void iigKernel(ars::Vec2d* mean1data, ars::Vec2d* mean2data, double sigma1, double sigma2, size_t kernelNum /*!!!! kernelNum = means.size() */, int fourierOrder, bool pnebiMode, ars::PnebiLUT& pnebiLUT, double* coefficients) {
    //    a.insertIsotropicGaussians(points, sigma);

    //TODO 1): 4 righe sotto: DA FARE NEL MAIN PRIMA DI CHIAMARE LA FUNZIONE; vengono fatte una tantum prima del for
    //    if (coeffs_.size() != 2 * arsfOrder_ + 2) {
    //        coeffs_.resize(2 * arsfOrder_ + 2);
    //    }
    //    std::fill(coeffs_.begin(), coeffs_.end(), 0.0);

    for (size_t i = 0; i < kernelNum; ++i) {
        for (size_t j = i + 1; j < kernelNum; ++j) {
            //            isotropicKer_.init(means[i], means[j], sigma);
            double dx, dy;

            ars::Vec2d vecI = mean1data[i];
            ars::Vec2d vecJ = mean2data[j];

            dx = vecJ.x - vecI.x;
            dy = vecJ.y - vecI.y;
            double phi = atan2(dy, dx);
            double sigmaValSq = sigma1 * sigma1 + sigma2 * sigma2;
            double lambdaSqNorm = 0.25 * (dx * dx + dy * dy) / sigmaValSq;


            //            isotropicKer_.updateFourier(arsfOrder_, coeffs_, w);
            double weight = 1.0 / (kernelNum * kernelNum);
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

            if (pnebiMode == false) {
                //                updateARSF2CoeffRecursDown(lambdaSqNorm, phi, w2, nFourier, coeffs);

                double cth2, sth2;
                cth2 = cos(2.0 * phi);
                sth2 = sin(2.0 * phi);
                //                updateARSF2CoeffRecursDown(lambda, cth2, sth2, factor, n, coeffs);


                //TODO 4): make pnebis a double*
                //can solve it with cuda/gpu malloc here? otherwise just pass the needed pointer to the function?
                //for now I just declare it here
                //                std::vector<double> pnebis(n + 1);
                double *pnebis;

                double sgn, cth, sth, ctmp, stmp;

                // Fourier Coefficients 
                //                if (coeffs.size() != 2 * n + 2) {
                //                    std::cerr << __FILE__ << "," << __LINE__ << ": invalid size of Fourier coefficients vector " << coeffs.size() << " should be " << (2 * n + 2) << std::endl;
                //                    return;
                //                }

                //                TODO 5): expand evaluatePnebiVector() below
                //                evaluatePnebiVector(n, lambda, pnebis);
                //ARS_PRINT(pnebis[0]);

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
            } else if (pnebiMode == true) {
                //                updateARSF2CoeffRecursDownLUT(lambdaSqNorm_, phi_, w2, nFourier, pnebiLut_, coeffs);
                double cth2, sth2;
                //fastCosSin(2.0 * phi, cth2, sth2); 
                cth2 = cos(2.0 * phi);
                sth2 = sin(2.0 * phi);

                //TODO 6): find a workaround for this pnebis vector; for now I just initialize here a double* pnebis;
                //                std::vector<double> pnebis(fourierOrder + 1); //prima riga della funzione omonima chiamata da dentro l'inline
                double* pnebis;
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
                //                {
                //                    double val, w;
                //
                //                    // Checkes LUT initialization
                //                    if (lut_.empty()) {
                //                        std::cerr << __FILE__ << "," << __LINE__ << ": empty LUT" << std::endl;
                //                        return;
                //                    }
                //
                //                    // Checks that 1) argument x is positive; 2) order of PNEBI has been computed;
                //                    if (x < 0.0) x = -x;
                //
                //                    // Checkes the size of vector (and adapt it if needed)
                //                    if (y.size() < orderMax_ + 1) {
                //                        y.resize(orderMax_ + 1);
                //                    }
                //
                //                    PnebiPoint tmp(0, x);
                //                    auto upper = std::upper_bound(lut_.begin(), lut_.end(), tmp);
                //                    auto lower = upper;
                //                    std::advance(lower, -1);
                //
                //                    if (upper == lut_.end()) {
                //                        std::copy(lower->y.begin(), lower->y.end(), y.begin());
                //                        //evaluatePnebiVector(orderMax_,x,y);
                //                    } else {
                //                        w = (x - lower->x) / (upper->x - lower->x);
                //                        for (int i = 0; i < lower->y.size() && i < upper->y.size() && i < y.size(); ++i) {
                //                            y[i] = (1.0 - w) * lower->y[i] + w * upper->y[i];
                //                        }
                //                    }
                //                }
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
}

int main(void) {
    ars::AngularRadonSpectrum2d ars1;
    ars::AngularRadonSpectrum2d ars2;
    thrust::device_vector<ars::Vec2d> acesPoints1;
    std::chrono::system_clock::time_point timeStart, timeStop;
    double sigma = 0.05;
    int fourierOrder = 20;

    ars1.setARSFOrder(fourierOrder);
    ars2.setARSFOrder(fourierOrder);

    rangeToPoint(acesRanges1, 180, -0.5 * M_PI, M_PI / 180.0 * 1.0, acesPoints1);
    //        acesPoints1.push_back(ars::Vector2::Zero());
    ars::Vec2d firstElement;
    //    firstElement.resetToZero(); 
    firstElement.x = 0.0;
    firstElement.y = 0.0;
    acesPoints1.push_back(firstElement);
    std::cout << "Number of input points: " << acesPoints1.size() << std::endl;
    //    for (int i = 0; i < acesPoints1.size(); ++i) {
    //        std::cout << i << "\t" << acesPoints1[i].x() << "\t" << acesPoints1[i].y() << std::endl;
    //    }

    //    const int trialNum = 1;
    ars1.initLUT(0.0001);
    ars1.setComputeMode(ars::ArsKernelIsotropic2d::ComputeMode::PNEBI_DOWNWARD);


    timeStart = std::chrono::system_clock::now();

    //kernel call
    //    ars1.insertIsotropicGaussians(acesPoints1, sigma);
    ars::Vec2d* kernelInput1 = thrust::raw_pointer_cast(acesPoints1.data());
    size_t kernelNum = acesPoints1.size(); //numero di punti in input
    bool pnebiMode = false;
    double *coefficientsArs1;
    ars::PnebiLUT pnebiLUT1; //LUT setup
    double lutPrecision = 0.001;
    pnebiLUT1.init(fourierOrder, lutPrecision); //LUT setup
    iigKernel << <1, 1 >> >(kernelInput1, kernelInput1, sigma, sigma, kernelNum, fourierOrder, pnebiMode, pnebiLUT1, coefficientsArs1);
    //end of kernel call

    timeStop = std::chrono::system_clock::now();
    double timeArs1 = (double) std::chrono::duration_cast<std::chrono::milliseconds>(timeStop - timeStart).count();
    cudaDeviceSynchronize();
    std::cout << "insertIsotropicGaussians() " << timeArs1 << " ms" << std::endl;


    std::cout << "\n------\n" << std::endl;

    ars2.setComputeMode(ars::ArsKernelIsotropic2d::ComputeMode::PNEBI_LUT);



    timeStart = std::chrono::system_clock::now();

    //kernel call
    //    ars2.insertIsotropicGaussians(acesPoints1, sigma);
    ars::Vec2d* kernelInput2 = thrust::raw_pointer_cast(acesPoints1.data());
    kernelNum = acesPoints1.size(); //for this dummy example atleast
    pnebiMode = true;
    double *coefficientsArs2;
    //    iigKernel << <1, 1 >> >(kernelInput1, kernelInput1, sigma, sigma, kernelNum, fourierOrder, pnebiMode, coefficientsArs2);
    //end of kernel call


    timeStop = std::chrono::system_clock::now();
    double timeArs2 = (double) std::chrono::duration_cast<std::chrono::milliseconds>(timeStop - timeStart).count();
    std::cout << "insertIsotropicGaussians() " << timeArs2 << " ms" << std::endl;



    std::cout << "\nARS Coefficients:\n";
    std::cout << "\ti \tDownward \tLUT\n";
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


    return 0;
}

void rangeToPoint(double* ranges, int num, double angleMin, double angleRes, thrust::device_vector<ars::Vec2d>& points) {
    ars::Vec2d p;
    for (int i = 0; i < num; ++i) {
        double a = angleMin + angleRes * i;
        ranges[i] * cos(a), ranges[i] * sin(a);
        points.push_back(p);
    }
}



