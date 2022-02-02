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
#include <unsupported/Eigen/FFT>
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
    	static constexpr double INV_SQRT_2_PI = 0.3989422804;   // 1 / sqrt(2.0 * M_PI)

        /**
         * Creates a flat kernel
         */
        ArsKernelAnisotropic2d();

        /**
         * Creates a kernel with the given number of Fourier coefficients.
         */
        ArsKernelAnisotropic2d(int nFourier);

        /**
         * Constructor of non-isotropic Kernel of Angular Randon Specturm associated to two Gaussian distributions. 
         * @param mean1 mean value of first gaussian
         * @param covar1 covariance matrix of first gaussian
         * @param mean2 mean value of second gaussian
         * @param covar2 covariance matrix of second gaussian
         * @param nFourier Fourier order to be computed
         */
        ArsKernelAnisotropic2d(const Vec2d& mean1, const Mat2d& covar1, const Vec2d& mean2, const Mat2d& covar2, int nFourier = 256);

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
        void init(const Vec2d& mean1, const Mat2d& covar1, const Vec2d& mean2, const Mat2d& covar2);
        
        /**
         * Sets the number of discrete intervals used in numerical integration of over 
         * the period M_PI.
         * @param nRes number of interval 
         */
//        void setIntervalNum(int nRes) {
//            nRes_ = nRes;
//        }

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

        void setFourierOrder(int nFourier) {
        	if (nFourier_ != nFourier) {
        		nFourier_ = nFourier;
        		initCosSinLut();
        	}
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
        void computeFourier(std::vector<double>& coeffs);
        
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
        //void updateFourier(std::vector<double>& coeffs);

    private:
        struct DataLut {
        	double varCos;
        	double varSin;
        	double meanConst;
        	double meanCos;
        	double meanSin;
        	std::vector<double> cosTh;
        	std::vector<double> sinTh;
        };

        int nFourier_;
        double muMod_;
        double muAng_;
        double sigmaMod_;
        double sigmaAng_;
        double sigmaDif_;
        Eigen::FFT<double> fft_;
        std::vector<double> kernelVal_;
        std::vector<std::complex<double> > freqvec_;
        DataLut lut_;

        void initCosSinLut();
    };

}
