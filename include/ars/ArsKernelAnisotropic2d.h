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
#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <ars/definitions.h>

namespace ars {

    /**
     * Class AnisotropicKernel computes and handles the kernel of Angular Radon Spectrum (ARS)
     * associated to a pair of Gaussian distributions in a Gaussian mixture. 
     * A Gaussian distribution is isotropic if its covariance matrix C is
     *      C = sigma^2 * I
     * i.e. C is diagonal with equal values of variance in all the directions. 
     * The ARS kernal of isotropic distributions can be expanded in Fourier series 
     * whose coefficients have a closed-form formula. 
     * The present anisotropic case has no known closed-form formula and, thus, 
     * the coeffients are computed through numerical integration. 
     * (Hence, the need to provide the 
     */
    class ArsKernelAnisotropic2d {
    public:
        /**
         * Creates a flat kernel
         */
        ArsKernelAnisotropic2d();

        /**
         * Constructor of non-isotropic Kernel of Angular Randon Specturm associated to two Gaussian distributions. 
         * @param mean1 mean value of first gaussian
         * @param covar1 covariance matrix of first gaussian
         * @param mean2 mean value of second gaussian
         * @param covar2 covariance matrix of second gaussian
         */
        ArsKernelAnisotropic2d(const Vector2& mean1, const Matrix2& covar1, const Vector2& mean2, const Matrix2& covar2);

        /**
         * Destructor. 
         */
        ~ArsKernelAnisotropic2d();

        /**
         * Computes the kernel parameters and save them. 
         * @param mean1 mean value of first gaussian
         * @param covar1 covariance matrix of first gaussian
         * @param mean2 mean value of second gaussian
         * @param covar2 covariance matrix of second gaussian
         */
        void init(const Vector2& mean1, const Matrix2& covar1, const Vector2& mean2, const Matrix2& covar2);
        
        /**
         * Sets the number of discrete intervals used in numerical integration of over 
         * the period M_PI.
         * @param nRes number of interval 
         */
        void setIntervalNum(int nRes) {
            nRes_ = nRes_;
        }

        /**
         * Returns the module of the sinusoidal numerator. 
         */
        double getMuModule() const {
            return muMod_;
        }

        /**
         * Return the phase of the sinusoidal numerator. 
         * @return 
         */
        double getMuPhase() const {
            return muAng_;
        }

        /**
         * Return the average value of the eigenvalues of covariance matrix. 
         * It also corresponds to the module of sinusoidal denominator.
         */
        double getVarianceModule() const {
            return sigmaMod_;
        }
        
        /**
         * Return the phase of the sinusoidal denominator.
         */
        double getVariancePhase() const {
            return sigmaAng_;
        }
        
        /**
         * Returns the "eccentricity" of the eigenvalues of covariance matrix. 
         */
        double getVariancePerc() const {
            return sigmaDif_;
        }

        /**
         * Copmputes the value of the non-isotropic kernel at a given angle. 
         * @param t
         * @return 
         */
        inline double value(double t) const {
            double var = sigmaMod_ * (1.0 + sigmaDif_ * cos(2.0 * t - 2.0 * sigmaAng_));
            double mean = 0.5 * muMod_ * muMod_ * (1.0 + cos(2.0 * t - 2.0 * muAng_));
            return exp(-0.5 * mean / var) / sqrt(2 * M_PI * var);
        }

        /**
         * Computes the coefficients of Fourier series expansion of the kernel. 
         * The series is truncated to n-th order. 
         * The M_PI period of the kernel is divided into k intervals to compute 
         * the coefficients using numeric integration. 
         * @param nFourier maximum order of the Fourier series 
         * @param nRes number of intervals used in numeric integration
         * @param coeffs the computed coefficients 
         */
        void computeFourier(int nFourier, std::vector<double>& coeffs) const;
        
        /**
         * Computes the coefficients of Fourier series expansion of the kernel 
         * and ADD them to the passed ones. (It does not reset them to zero!)
         * The series is truncated to n-th order. 
         * The M_PI period of the kernel is divided into k intervals to compute 
         * the coefficients using numeric integration. 
         * @param nFourier maximum order of the Fourier series 
         * @param nRes number of intervals used in numeric integration
         * @param coeffs the computed coefficients 
         */
        void updateFourier(int nFourier, std::vector<double>& coeffs) const;

    private:
        double muMod_;
        double muAng_;
        double sigmaMod_;
        double sigmaAng_;
        double sigmaDif_;
        int nRes_;
    };

}
