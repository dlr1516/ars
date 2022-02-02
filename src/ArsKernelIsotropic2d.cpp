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
#include <ars/ArsKernelIsotropic2d.h>

namespace ars {

    ArsKernelIsotropic2d::ArsKernelIsotropic2d() : lambdaSqNorm_(0.0), sigmaValSq_(1.0), phi_(0.0), pnebiLut_(40, 0.001), mode_(ComputeMode::PNEBI_LUT) {
    }

    ArsKernelIsotropic2d::ArsKernelIsotropic2d(const Vec2d& mean1, const Vec2d& mean2, double sigma)
    : pnebiLut_(40, 0.001), mode_(ComputeMode::PNEBI_LUT) {
        init(mean1, mean2, sigma);
    }

    ArsKernelIsotropic2d::ArsKernelIsotropic2d(const Vec2d& mean1, const Vec2d& mean2, double sigma1, double sigma2)
    : pnebiLut_(40, 0.001), mode_(ComputeMode::PNEBI_LUT) {
        init(mean1, mean2, sigma1, sigma2);
    }

    ArsKernelIsotropic2d::~ArsKernelIsotropic2d() {
    }

    void ArsKernelIsotropic2d::init(const Vec2d& mean1, const Vec2d& mean2, double sigma) {
        init(mean1, mean2, sigma, sigma);
    }

    void ArsKernelIsotropic2d::init(const Vec2d& mean1, const Vec2d& mean2, double sigma1, double sigma2) {
        double dx, dy;

        dx = mean2.data_[0] - mean1.data_[0];
        dy = mean2.data_[1] - mean1.data_[1];
        phi_ = atan2(dy, dx);
        sigmaValSq_ = sigma1 * sigma1 + sigma2 * sigma2;
        lambdaSqNorm_ = 0.25 * (dx * dx + dy * dy) / sigmaValSq_;
    }

    void ArsKernelIsotropic2d::initPnebiLut(int n, double tol) {
        pnebiLut_.init(n, tol);
    }

    void ArsKernelIsotropic2d::computeFourier(int nFourier, std::vector<double>& coeffs) {
        if (coeffs.size() != 2 * nFourier + 2) {
            coeffs.resize(2 * nFourier + 2);
        }
        std::fill(coeffs.begin(), coeffs.end(), 0.0);

        updateFourier(nFourier, coeffs);
    }
    
    void ArsKernelIsotropic2d::updateFourier(int nFourier, std::vector<double>& coeffs) {
        updateFourier(nFourier, coeffs, 1.0);
    }

    void ArsKernelIsotropic2d::updateFourier(int nFourier, std::vector<double>& coeffs, double weight) {
        double w = weight / sqrt(2.0 * M_PI * sigmaValSq_);
        
        if (coeffs.size() != 2 * nFourier + 2) {
            coeffs.resize(2 * nFourier + 2);
        }
        
        if (pnebiLut_.getOrderMax() < nFourier) {
            ARS_ERROR("LUT not initialized to right order. Initialized now.");
            pnebiLut_.init(nFourier, 0.0001);
        }

        if (mode_ == ComputeMode::PNEBI_DOWNWARD) {
            //updateARSF2CoeffRecursDown(lambda, ux * ux - uy*uy, 2.0 * ux * uy, 1.0, arsfOrder_, coeffs_);
            updateARSF2CoeffRecursDown(lambdaSqNorm_, phi_, w, nFourier, coeffs);
        } else if (mode_ == ComputeMode::PNEBI_LUT) {
            //updateARSF2CoeffRecursDownLUT(lambda, ux * ux - uy*uy, 2.0 * ux * uy, 1.0, arsfOrder_, pnebiLut_, coeffs_);
            updateARSF2CoeffRecursDownLUT(lambdaSqNorm_, phi_, w, nFourier, pnebiLut_, coeffs);
        }
    }

    // ----------------------------------------------------
    // PRIVATE MEMBERS
    // ----------------------------------------------------

    std::array<std::string, 2> const ArsKernelIsotropic2d::MODE_NAME{"PNEBI_DOWNWARD", "PNEBI_LUT"};


} // end of namespace

