/**
 * ARS - Angular Radon Spectrum 
 * Copyright (C) 2020 Dario Lodi Rizzini.
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
#include <ars/ArsKernelAnisotropic2d.h>
#include <ars/utils.h>

#include <ars/functions.h>

namespace ars {

    ArsKernelAnisotropic2d::ArsKernelAnisotropic2d()
    : muMod_(0.0), muAng_(0.0), sigmaMod_(0.0), sigmaAng_(0.0), sigmaDif_(0.0), fft_(), freqvec_() {
    }

    ArsKernelAnisotropic2d::ArsKernelAnisotropic2d(const Vector2& mean1, const Matrix2& covar1, const Vector2& mean2, const Matrix2& covar2)
    : fft_(), freqvec_() {
        init(mean1, covar1, mean2, covar2);
    }

    ArsKernelAnisotropic2d::~ArsKernelAnisotropic2d() {
    }

    void ArsKernelAnisotropic2d::init(const Vector2& mean1, const Matrix2& covar1, const Vector2& mean2, const Matrix2& covar2) {
        Vector2 mu12;
        Matrix2 sigma12;
        double a, b, lmax, lmin, c, s;

        mu12 = mean2 - mean1;
        muMod_ = mu12.norm();
        muAng_ = atan2(mu12(1), mu12(0));
        sigma12 = covar1 + covar2;

        // Diagonalizes sigma12
        diagonalize(sigma12, lmin, lmax, sigmaAng_);

        //        a = 0.5 * (sigma12(1, 1) - sigma12(0, 0));
        //        b = 0.5 * (sigma12(0, 1) + sigma12(1, 0));
        //        ARS_VARIABLE2(a, b);
        //
        //        sigmaAng_ = 0.5 * atan2(-b, a);
        //
        //        c = cos(sigmaAng_);
        //        s = sin(sigmaAng_);
        //        lmax = sigma12(0, 0) * c * c + sigma12(1, 1) * s * s + (sigma12(0, 1) + sigma12(1, 0)) * c * s;
        //        lmin = sigma12(0, 0) * s * s + sigma12(1, 1) * c * c - (sigma12(0, 1) + sigma12(1, 0)) * c * s;
        //        ARS_VARIABLE3(sigmaAng_, lmax, lmin);
        //
        //        if (lmax < lmin) {
        //            sigmaAng_ += 0.5 * M_PI;
        //            std::swap(lmax, lmin);
        //            ARS_PRINT("lmin " << lmin << " < lmax " << lmax << ": swap, sigmaAng_ + PI/2: " << sigmaAng_);
        //        }

        sigmaMod_ = 0.5 * (lmax + lmin);
        sigmaDif_ = (lmax - lmin) / (lmax + lmin);

        //        ARS_PRINT("muMod_ " << muMod_ << ", muAng_[rad] " << muAng_ << " [deg] " << (180.0 / M_PI * muAng_) << "\n"
        //                << "sigmaMod_ " << sigmaMod_ << ", sigmaAng_[rad] " << sigmaAng_ << " [deg] " << (180.0 / M_PI * sigmaAng_)
        //                << ", sigmaDif_ " << sigmaDif_ << "\n");
    }

    void ArsKernelAnisotropic2d::computeFourier(int nFourier, std::vector<double>& coeffs) {
        if (coeffs.size() != 2 * nFourier + 2) {
            coeffs.resize(2 * nFourier + 2);
        }
        std::fill(coeffs.begin(), coeffs.end(), 0.0);

        updateFourier(nFourier, coeffs);
    }

    void ArsKernelAnisotropic2d::updateFourier(int nFourier, std::vector<double>& coeffs) {
        std::vector<double> kernelVal(nFourier);
        double sumCos, sumSin, cosCurr, cosNext, cosIncr, sinCurr, sinNext, sinIncr;
        double dt, factor;
        //double h = 0.5 * dt / M_PI;

        ARS_ASSERT(nFourier > 0);

        // Evaluates the kernel function at given intervals
        dt = M_PI / nFourier;
        for (int i = 0; i < nFourier; ++i) {
            kernelVal[i] = value(dt * i);
        }
        // Alternative and possibly more efficient evaluation
//        double varCos = sigmaMod_ * sigmaDif_ * cos(2.0 * sigmaAng_);
//        double varSin = sigmaMod_ * sigmaDif_ * sin(2.0 * sigmaAng_);
//        double meanConst = 0.5 * muMod_ * muMod_;
//        double meanCos = meanConst * cos(2.0 * muAng_);
//        double meanSin = meanConst * sin(2.0 * muAng_);
//        double varInv, cosTh, sinTh;
//        for (int i = 0; i < nFourier; ++i) {
//        	//cosTh = cos(dt * i);
//        	//sinTh = sin(dt * i);
//        	ars::fastCosSin(dt*i, cosTh, sinTh);
//        	varInv = 1.0 / (sigmaMod_ + varCos * cosTh + varSin * sinTh);
//        	kernelVal[i] = INV_SQRT_2_PI * sqrt(varInv) * exp(-0.5 * (meanConst + meanCos * cosTh + meanSin * sinTh) * varInv);
//        	ARS_ASSERT(fabs(kernelVal[i] - value(dt * i)) < 1e-6)
//        }


        // Evaluates each of the Fourier coefficients
        if (coeffs.size() != 2 * nFourier + 2) {
            coeffs.resize(2 * nFourier + 2);
        }

        // Computation of the 0-order coefficient (only for cosine part) using 
        // numerical integration (trapezoid approximation)
        //        for (int i = 0; i < nRes_; ++i) {
        //            coeffs[0] += (kernelVal[i] + kernelVal[i + 1]) * h;
        //        }

        // Computes the Fourier coefficients of orders greater or equal to 1 
        // according to trapezoidal rule integration
        //        cosIncr = cos(2.0 * dt);
        //        sinIncr = sin(2.0 * dt);
        //        h = 2.0 * h;
        //        for (int k = 1; k <= nFourier; ++k) {
        //            sumCos = 0.0;
        //            sumSin = 0.0;
        //            cosCurr = 1.0;
        //            sinCurr = 0.0;
        //            for (int i = 0; i < nRes_; ++i) {
        //                fastCosSin(2.0 * k * dt * (i + 1), cosNext, sinNext);
        //                //                cosNext = cos(2.0 * k * dt * (i + 1));
        //                //                sinNext = sin(2.0 * k * dt * (i + 1));
        //                //                cosNext = cosCurr * cosIncr - sinCurr * sinIncr;
        //                //                sinNext = cosCurr * sinIncr + sinCurr * cosIncr;
        //                sumCos += (kernelVal[i] * cosCurr + kernelVal[i + 1] * cosNext) * h;
        //                sumSin += (kernelVal[i] * sinCurr + kernelVal[i + 1] * sinNext) * h;
        //                cosCurr = cosNext;
        //                sinCurr = sinNext;
        //            }
        //            coeffs[2 * k] = sumCos;
        //            coeffs[2 * k + 1] = sumSin;
        //        }
        //        
        //        std::vector<double> coeffsFft;
        //        fft(kernelVal, coeffsFft, nFourier);

        //        ARS_PRINT("Compare integral Fourier and FFT: coeffs.size() " << coeffs.size() << ", coeffsFft.size() " << coeffsFft.size());
        //        for (int i = 0; i < coeffs.size() && i < coeffsFft.size(); ++i) {
        //            std::cout << "  i " << i << ": \t" << coeffs[i] << " \t" << coeffsFft[i] << "\n";
        //        }
        fft(kernelVal, coeffs, nFourier);

//		fft_.fwd(freqvec_, kernelVal);

		//ARS_PRINT("fft: input size " << funIn.size() << ", output complex size " << freqvec.size());

//		coeffs[0] = freqvec_[0].real() / nFourier;
//		coeffs[1] = 0.0;
//		factor = 2.0 / nFourier;
//		for (int i = 1; i <= nFourier; ++i) {
//			coeffs[2 * i] = factor * freqvec_[i].real();
//			coeffs[2 * i + 1] = -factor * freqvec_[i].imag();
//		}
    }

}
