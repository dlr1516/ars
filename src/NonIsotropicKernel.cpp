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
#include <ars/NonIsotropicKernel.h>

namespace ars {

    NonIsotropicKernel::NonIsotropicKernel()
    : muMod_(0.0), muAng_(0.0), sigmaMod_(0.0), sigmaAng_(0.0), sigmaDif_(0.0) {
    }

    NonIsotropicKernel::NonIsotropicKernel(const Vector2& mean1, const Matrix2& covar1, const Vector2& mean2, const Matrix2& covar2) {
        init(mean1, covar1, mean2, covar2);
    }

    NonIsotropicKernel::~NonIsotropicKernel() {
    }

    void NonIsotropicKernel::init(const Vector2& mean1, const Matrix2& covar1, const Vector2& mean2, const Matrix2& covar2) {
        Vector2 mu12;
        Matrix2 sigma12;
        double a, b, l1, l2, c, s;

        mu12 = mean2 - mean1;
        muMod_ = mu12.norm();
        muAng_ = atan2(mu12(1), mu12(0));
        sigma12 = covar1 + covar2;

        // Diagonalizes sigma12
        a = 0.5 * (sigma12(1, 1) - sigma12(0, 0));
        b = 0.5 * (sigma12(0, 1) + sigma12(1, 0));
        ARS_VARIABLE2(a, b);

        sigmaAng_ = 0.5 * atan2(-b, a);

        c = cos(sigmaAng_);
        s = sin(sigmaAng_);
        l1 = sigma12(0, 0) * c * c + sigma12(1, 1) * s * s + (sigma12(0, 1) + sigma12(1, 0)) * c * s;
        l2 = sigma12(0, 0) * s * s + sigma12(1, 1) * c * c - (sigma12(0, 1) + sigma12(1, 0)) * c * s;
        ARS_VARIABLE3(sigmaAng_, l1, l2);

        if (l1 < l2) {
            sigmaAng_ += 0.5 * M_PI;
            std::swap(l1, l2);
            ARS_PRINT("l1 " << l1 << " < l2 " << l2 << ": swap, sigmaAng_ + PI/2: " << sigmaAng_);
        }

        sigmaMod_ = 0.5 * (l1 + l2);
        sigmaDif_ = (l1 - l2) / (l1 + l2);

        ARS_PRINT("muMod_ " << muMod_ << ", muAng_[rad] " << muAng_ << " [deg] " << (180.0 / M_PI * muAng_) << "\n"
                << "sigmaMod_ " << sigmaMod_ << ", sigmaAng_[rad] " << sigmaAng_ << " [deg] " << (180.0 / M_PI * sigmaAng_)
                << ", sigmaDif_ " << sigmaDif_ << "\n");
    }

    void NonIsotropicKernel::computeFourier(int nFourier, int nRes, std::vector<double>& coeffs) const {
        std::vector<double> kernelVal(nRes + 1);
        double sumCos, sumSin, cosCurr, cosNext, cosIncr, sinCurr, sinNext, sinIncr;
        double dt = M_PI / nRes;
        double dt2 = 0.5 * dt;

        // Evaluates the kernel function at given intervals
        for (int i = 0; i <= nRes; ++i) {
            kernelVal[i] = value(dt * i);
        }

        // Evaluates each of the Fourier coefficients
        coeffs.resize(2 * nFourier + 2, 0.0);

        // Computation of the 0-order coefficient (only for cosine part) using 
        // numerical integration (trapezoid approximation)
        for (int i = 0; i < nRes; ++i) {
            coeffs[0] += (kernelVal[i] + kernelVal[i + 1]) * dt2;
        }

        // Computes the Fourier coefficients of orders greater or equal to 1 
        // according to trapezoidal rule integration
        cosIncr = cos(2.0 * dt);
        sinIncr = sin(2.0 * dt);
        for (int k = 1; k <= nFourier; ++k) {
            sumCos = 0.0;
            sumSin = 0.0;
            cosCurr = 1.0;
            sinCurr = 0.0;
            for (int i = 0; i < nRes; ++i) {
                //cosNext = cos(2.0 * dt * (i+1));
                //sinNext = sin(2.0 * dt * (i+1));
                cosNext = cosCurr * cosIncr - sinCurr * sinIncr;
                sinNext = cosCurr * sinIncr + sinCurr * cosIncr;
                sumCos += (kernelVal[i] * cosCurr + kernelVal[i + 1] * cosNext) * dt2;
                sumSin += (kernelVal[i] * sinCurr + kernelVal[i + 1] * sinNext) * dt2;
                cosCurr = cosNext;
                sinCurr = sinNext;
            }
            coeffs[2 * k] = sumCos;
            coeffs[2 * k + 1] = sumSin;
        }
    }

}
