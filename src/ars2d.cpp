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
#include <ars/ars2d.h>
#include <boost/numeric/interval.hpp>
#include <queue>

//#include <omp.h>

namespace ars {

    // --------------------------------------------------------
    // ARSF FUNCTIONS
    // --------------------------------------------------------

    /** Computes the contribution of term (lambda,phi) to the Fourier coefficients of ARS
     * computed with downward recursion of modified Bessel functions.
     */
    void updateARSF2CoeffRecursDown(double lambda, double cth2, double sth2, double factor, int n, std::vector<double>& coeffs) {
        std::vector<double> pnebis(n + 1);
        double sgn, cth, sth, ctmp, stmp;

        // Fourier Coefficients 
        if (coeffs.size() != 2 * n + 2) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid size of Fourier coefficients vector " << coeffs.size() << " should be " << (2 * n + 2) << std::endl;
            return;
        }

        evaluatePnebiVector(n, lambda, pnebis);
        //ARS_PRINT(pnebis[0]);
        coeffs[0] += 0.5 * factor * pnebis[0];
        sgn = -1.0;
        cth = cth2;
        sth = sth2;
        for (int k = 1; k <= n; ++k) {
            coeffs[2 * k] += factor * pnebis[k] * sgn * cth;
            coeffs[2 * k + 1] += factor * pnebis[k] * sgn * sth;
            sgn = -sgn;
            ctmp = cth2 * cth - sth2 * sth;
            stmp = sth2 * cth + cth2 * sth;
            cth = ctmp;
            sth = stmp;
        }
    }

    void updateARSF2CoeffRecursDownLUT(double lambda, double cth2, double sth2, double factor, int n, const PnebiLUT& pnebiLUT, std::vector<double>& coeffs) {
        std::vector<double> pnebis(n + 1);
        double sgn, cth, sth, ctmp, stmp;

        // Fourier Coefficients 
        if (coeffs.size() != 2 * n + 2 || pnebiLUT.getOrderMax() < n) {
            std::cerr << __FILE__ << "," << __LINE__ << ": one of these conditions failed:"
                    << "\n  size of Fourier coefficients vector " << coeffs.size() << " should be " << (2 * n + 2)
                    << "\n  LUT max order is " << pnebiLUT.getOrderMax() << " >= " << n
                    << std::endl;
            return;
        }

        pnebiLUT.eval(lambda, pnebis);
        //ARS_PRINT(pnebis[0]);
        coeffs[0] += 0.5 * factor * pnebis[0];
        sgn = -1.0;
        cth = cth2;
        sth = sth2;
        for (int k = 1; k <= n; ++k) {
            coeffs[2 * k] += pnebis[k] * factor * sgn * cth;
            coeffs[2 * k + 1] += pnebis[k] * factor * sgn * sth;
            sgn = -sgn;
            ctmp = cth2 * cth - sth2 * sth;
            stmp = sth2 * cth + cth2 * sth;
            cth = ctmp;
            sth = stmp;
        }
    }

    void computeFourierCorr(const std::vector<double>& fourierSrc, const std::vector<double>& fourierDst, std::vector<double>& fourierCor) {
        int n;

        if (fourierSrc.size() % 2 != 0 || fourierDst.size() % 2 != 0 || fourierSrc.size() != fourierDst.size()) {
            std::cerr << __FILE__ << "," << __LINE__ << ": the number of ARSF coefficients must be even and equal: fourierSrc " << fourierSrc.size()
                    << ", fourierDst " << fourierDst.size() << std::endl;
            return;
        }
        n = (fourierSrc.size() / 2) - 1;

        // Resizes the correlation coefficients, if required
        if (fourierCor.size() != 2 * n + 2) {
            fourierCor.resize(2 * n + 2);
        }

        // Computes the coefficients
        for (int k = 0; k <= n; ++k) {
            fourierCor[2 * k] = 0.5 * (fourierSrc[2 * k] * fourierDst[2 * k] + fourierSrc[2 * k + 1] * fourierDst[2 * k + 1]);
            fourierCor[2 * k + 1] = 0.5 * (fourierSrc[2 * k] * fourierDst[2 * k + 1] - fourierSrc[2 * k + 1] * fourierDst[2 * k]);
        }
    }

    // --------------------------------------------------------
    // ARS 2D CLASS
    // --------------------------------------------------------

    std::array<std::string, 2> const AngularRadonSpectrum2d::MODE_NAME{"PNEBI_DOWNWARD", "PNEBI_LUT"};

    AngularRadonSpectrum2d::AngularRadonSpectrum2d()
    : coeffs_(), arsfOrder_(0),
    thetaToll_(M_PI / 180.0 * 0.5), threadNumOMP_(4), pnebiLut_(), mode_(PNEBI_LUT) {
    }

    AngularRadonSpectrum2d::AngularRadonSpectrum2d(const std::vector<double>& coeffs)
    : coeffs_(coeffs), arsfOrder_(0),
    thetaToll_(M_PI / 180.0 * 0.5), threadNumOMP_(4), pnebiLut_(), mode_(PNEBI_LUT) {
    }

    AngularRadonSpectrum2d::~AngularRadonSpectrum2d() {
    }

    void AngularRadonSpectrum2d::insertIsotropicGaussians(const Vector2Vector& means, double sigma) {
        int kernelNum = means.size();
        //std::cout << "kernelNum " << kernelNum << ", mode_ " << mode_ << " " << MODE_NAME[mode_] << std::endl;

        if (pnebiLut_.getOrderMax() < arsfOrder_) {
            std::cerr << __FILE__ << "," << __LINE__ << ": LUT not initialized to right order. Initialized now." << std::endl;
            pnebiLut_.init(arsfOrder_, 0.0001);
        }

        std::fill(coeffs_.begin(), coeffs_.end(), 0.0);
//#pragma omp parallel num_threads(threadNumOMP_) shared(means,sigmas,kernelNum) 
        double dx, dy, sigma2, lambda, phi, scale, ux, uy;
        for (int i = 0; i < kernelNum; ++i) {
            for (int j = i + 1; j < kernelNum; ++j) {
                dx = means[i].x() - means[j].x();
                dy = means[i].y() - means[j].y();
                sigma2 = 2.0 * sigma * sigma;
                lambda = (dx * dx + dy * dy);
                phi = atan2(dy, dx);
                scale = 1.0 / sqrt(lambda);
                ux = dx * scale;
                uy = dy * scale;
                lambda = lambda / (2.0 * sigma2);
//#pragma omp atomic
                //std::cout << "i " << i << ", j " << j << ": lambda " << lambda << ", phi " << phi << std::endl;
                if (mode_ == PNEBI_DOWNWARD) {
                    //updateARSF2CoeffRecursDown(lambda, ux * ux - uy*uy, 2.0 * ux * uy, 1.0, arsfOrder_, coeffs_);
                    updateARSF2CoeffRecursDown(lambda, phi, 1.0, arsfOrder_, coeffs_);
                } else if (mode_ == PNEBI_LUT) {
                    //updateARSF2CoeffRecursDownLUT(lambda, ux * ux - uy*uy, 2.0 * ux * uy, 1.0, arsfOrder_, pnebiLut_, coeffs_);
                    updateARSF2CoeffRecursDownLUT(lambda, phi, 1.0, arsfOrder_, pnebiLut_, coeffs_);
                }
            }
        }
    }

    void AngularRadonSpectrum2d::insertIsotropicGaussians(const Vector2Vector& means, const std::vector<double>& sigmas) {
        int kernelNum = means.size();

        if (kernelNum != sigmas.size()) {
            std::cerr << __FILE__ << "," << __LINE__ << ": inconsistent vector sizes: found " << means.size()
                    << " mean values and " << sigmas.size() << " standard deviations" << std::endl;
            return;
        }

        if (pnebiLut_.getOrderMax() < arsfOrder_) {
            //std::cerr << __FILE__ << "," << __LINE__ << ": LUT not initialized to right order. Initialized now." << std::endl;
            pnebiLut_.init(arsfOrder_, 0.005);
        }

        std::fill(coeffs_.begin(), coeffs_.end(), 0.0);
#pragma omp parallel num_threads(threadNumOMP_) shared(means,sigmas,kernelNum) 
        double dx, dy, sigma2, lambda, scale, ux, uy;
        for (int i = 0; i < kernelNum; ++i) {
            for (int j = i + 1; j < kernelNum; ++j) {
                dx = means[i].x() - means[j].x();
                dy = means[i].y() - means[j].y();
                sigma2 = sigmas[i] * sigmas[i] + sigmas[j] * sigmas[j];
                lambda = (dx * dx + dy * dy);
                scale = 1.0 / sqrt(lambda);
                ux = dx * scale;
                uy = dy * scale;
                lambda = lambda / (2.0 * sigma2);
#pragma omp atomic
                if (mode_ == PNEBI_DOWNWARD) {
                    updateARSF2CoeffRecursDown(lambda, ux * ux - uy*uy, 2.0 * ux * uy, 1.0, arsfOrder_, coeffs_);
                } else {
                    updateARSF2CoeffRecursDownLUT(lambda, ux * ux - uy*uy, 2.0 * ux * uy, 1.0, arsfOrder_, pnebiLut_, coeffs_);
                }
            }
        }
    }

    //    void AngularRadonSpectrum2d::initARSFRecursDown() {
    //        std::fill(coeffs_.begin(), coeffs_.end(), 0.0);
    //        //std::cout << __FILE__ << "," << __LINE__ << ": \n" << coeffsRecursDown_.transpose() << std::endl;
    //        //  std::cout << __FILE__ << "," << __LINE__ << ": coeffsRecursDown_\n" << std::endl;
    //        for (auto& p : pairData_) {
    //            updateARSF2CoeffRecursDown(p.lambda, p.phi, 1.0, arsfOrder_, coeffs_);
    //        }
    //        //        coeffsRecursDown_(0) += 0.5 * pointNum_;
    //        //        coeffsRecursDown_ = coeffsRecursDown_ / (sigma_ * sqrt(M_PI));
    //    }
    //
    //    void AngularRadonSpectrum2d::initARSFRecursDownLUT() {
    //        // Check PNEBI function LUT
    //        if (pnebiLut_.getOrderMax() < arsfOrder_) {
    //            std::cerr << __FILE__ << "," << __LINE__ << ": LUT not initialized to right order. Initialized now." << std::endl;
    //            pnebiLut_.init(arsfOrder_, 0.005);
    //        }
    //
    //        // Checking Fourier coefficients vector
    //        if (coeffs_.size() != 2 * arsfOrder_ + 2) {
    //            std::cerr << __FILE__ << "," << __LINE__ << ": size of Fourier coefficients vector " << coeffs_.size()
    //                    << " should be " << (2 * arsfOrder_ + 2) << ": resized" << std::endl;
    //            coeffs_.resize(2 * arsfOrder_ + 2);
    //        }
    //
    //        // Variables
    //        std::vector<double> pnebis(arsfOrder_ + 1);
    //        std::fill(coeffs_.begin(), coeffs_.end(), 0.0);
    //        for (auto& p : pairData_) {
    //            updateARSF2CoeffRecursDownLUT(p.lambda, p.phi, 1.0, arsfOrder_, pnebiLut_, coeffs_);
    //        }
    //    }

    /** Evaluates the ARSF using the coefficients obtained from downward recursion. 
     */
    double AngularRadonSpectrum2d::eval(double theta) const {
        evaluateFourier(coeffs_, 2.0 * theta);
    }


} // end of namespace

