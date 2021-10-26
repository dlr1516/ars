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

	ArsKernelAnisotropic2d::ArsKernelAnisotropic2d() :
			nFourier_(256), muMod_(0.0), muAng_(0.0), sigmaMod_(0.0), sigmaAng_(0.0), sigmaDif_(0.0), fft_(), kernelVal_(), freqvec_() {
		initCosSinLut();
	}

	ArsKernelAnisotropic2d::ArsKernelAnisotropic2d(int nFourier) :
			nFourier_(nFourier), muMod_(0.0), muAng_(0.0), sigmaMod_(0.0), sigmaAng_(0.0), sigmaDif_(0.0), fft_(), kernelVal_(), freqvec_() {
		initCosSinLut();
	}

	ArsKernelAnisotropic2d::ArsKernelAnisotropic2d(const Vector2 &mean1, const Matrix2 &covar1, const Vector2 &mean2, const Matrix2 &covar2, int nFourier) :
			nFourier_(nFourier), fft_(), kernelVal_(), freqvec_() {
		init(mean1, covar1, mean2, covar2);
		initCosSinLut();
	}

	ArsKernelAnisotropic2d::~ArsKernelAnisotropic2d() {
	}

	void ArsKernelAnisotropic2d::init(const Vector2 &mean1, const Matrix2 &covar1, const Vector2 &mean2, const Matrix2 &covar2) {
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

		// Initizalizes LUT
		lut_.varCos = sigmaMod_ * sigmaDif_ * cos(2.0 * sigmaAng_);
		lut_.varSin = sigmaMod_ * sigmaDif_ * sin(2.0 * sigmaAng_);
		lut_.meanConst = 0.5 * muMod_ * muMod_;
		lut_.meanCos = lut_.meanConst * cos(2.0 * muAng_);
		lut_.meanSin = lut_.meanConst * sin(2.0 * muAng_);
	}

//	void ArsKernelAnisotropic2d::computeFourier(std::vector<double> &coeffs) {
//		if (coeffs.size() != 2 * nFourier_ + 2) {
//			coeffs.resize(2 * nFourier_ + 2);
//		}
//		std::fill(coeffs.begin(), coeffs.end(), 0.0);
//
//		updateFourier(coeffs);
//	}

	void ArsKernelAnisotropic2d::computeFourier(std::vector<double> &coeffs) {
		//double sumCos, sumSin, cosCurr, cosNext, cosIncr, sinCurr, sinNext, sinIncr;
		double dt, factor, varInv;

		ARS_ASSERT(nFourier_ > 0);

		// Evaluates each of the Fourier coefficients
		if (coeffs.size() != 2 * nFourier_ + 2) {
			coeffs.resize(2 * nFourier_ + 2);
		}

		// Evaluates the kernel function at given intervals
//        dt = M_PI / nFourier_;
//        for (int i = 0; i < nFourier; ++i) {
//            kernelVal[i] = value(dt * i);
//        }
		// Alternative and faster evaluation using LUT
        //  var = sigmaMod_ * (1.0 + sigmaDif_ * cos(2.0 * t - 2.0 * sigmaAng_))
		//      = sigmaMod_ + sigmaMod_ * sigmaDif_ * (cos(2.0 * t) * cos(2.0 * sigmaAng_) + sin(2.0 * t) * sin(2.0 * sigmaAng_))
		//      = sigmaMod_ + (sigmaMod_ * sigmaDif_ * cos(2.0 * sigmaAng_)) * cos(2.0 * t) + (sigmaMod_ * sigmaDif_ * sin(2.0 * sigmaAng_)) * sin(2.0 * t)
		//      = sigmaMod_ + lut_.varCos * cosTh[i] + lut_.varSin * lut_.sinTh[i]
		for (int i = 0; i < nFourier_; ++i) {
			varInv = 1.0 / (sigmaMod_ + lut_.varCos * lut_.cosTh[i] + lut_.varSin * lut_.sinTh[i]);
			kernelVal_[i] = INV_SQRT_2_PI * sqrt(varInv) * exp(-0.5 * (lut_.meanConst + lut_.meanCos * lut_.cosTh[i] + lut_.meanSin * lut_.sinTh[i]) * varInv);
		}

		//fft(kernelVal, coeffs, nFourier);

		fft_.fwd(freqvec_, kernelVal_);

		coeffs[0] = freqvec_[0].real() / nFourier_;
		coeffs[1] = 0.0;
		factor = 2.0 / nFourier_;
		for (int i = 1; i < nFourier_; ++i) {
			coeffs[2 * i] = factor * freqvec_[i].real();
			coeffs[2 * i + 1] = -factor * freqvec_[i].imag();
		}
	}

	void ArsKernelAnisotropic2d::initCosSinLut() {
		double dt = 2.0 * M_PI / nFourier_;
		lut_.cosTh.resize(nFourier_);
		lut_.sinTh.resize(nFourier_);
		for (int i = 0; i < nFourier_; ++i) {
			lut_.cosTh[i] = cos(dt * i);
			lut_.sinTh[i] = sin(dt * i);
		}
		kernelVal_.resize(nFourier_);
	}

}
