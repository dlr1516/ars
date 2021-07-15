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
#ifndef ISOTROPICKERNEL_H
#define ISOTROPICKERNEL_H

#include <iostream>
#include <vector>
#include <cmath>
#include <type_traits>   // for std::underlying_type_t()
#include <ars/definitions.h>
#include <ars/functions.h>


namespace ars {

    // --------------------------------------------------------
    // ARS-2D FUNCTIONS FOR ISOTROPIC KERNEL
    // --------------------------------------------------------

    /** Computes the contribution of term (lambda,phi) to the Angular Radon Spectrum Fourier (ARSF) 
     * coefficients computed with downward recursion of modified Bessel functions. 
     * The vector of coefficients coeffs[i] are used in Fourier series:
     *   S(x) = \sum_{i=0}^{n} ( coeffs[2*i] * cos(2*i*x) + coeffs[2*i+1] * sin(2*i*x) )
     */
    void updateARSF2CoeffRecursDown(double lambda, double cth2, double sth2, double factor, int n, std::vector<double>& coeffs);

    /** Computes the contribution of term (lambda,phi) to the Angular Radon Spectrum Fourier (ARSF) 
     * coefficients computed with downward recursion of modified Bessel functions. 
     * The vector of coefficients coeffs[i] are used in Fourier series:
     *   S(x) = \sum_{i=0}^{n} ( coeffs[2*i] * cos(2*i*x) + coeffs[2*i+1] * sin(2*i*x) )
     */
    inline void updateARSF2CoeffRecursDown(double lambda, double phi, double factor, int n, std::vector<double>& coeffs) {
        double cth2, sth2;
        cth2 = cos(2.0 * phi);
        sth2 = sin(2.0 * phi);
        updateARSF2CoeffRecursDown(lambda, cth2, sth2, factor, n, coeffs);
    }

    /** Computes the contribution of term (lambda,phi) to the Angular Radon Spectrum Fourier (ARSF) 
     * coefficients computed with downward recursion of modified Bessel functions.
     * Values are computed using LUT!
     * The vector of coefficients coeffs[i] are used in Fourier series:
     *   S(x) = \sum_{i=0}^{n} ( coeffs[2*i] * cos(2*i*x) + coeffs[2*i+1] * sin(2*i*x) )
     */
    void updateARSF2CoeffRecursDownLUT(double lambda, double cth2, double sth2, double factor, int n, const PnebiLUT& pnebiLUT, std::vector<double>& coeffs);

    /** Computes the contribution of term (lambda,phi) to the Angular Radon Spectrum Fourier (ARSF) 
     * coefficients computed with downward recursion of modified Bessel functions.
     * Values are computed using LUT!
     * The vector of coefficients coeffs[i] are used in Fourier series:
     *   S(x) = \sum_{i=0}^{n} ( coeffs[2*i] * cos(2*i*x) + coeffs[2*i+1] * sin(2*i*x) )
     */
    inline void updateARSF2CoeffRecursDownLUT(double lambda, double phi, double factor, int n, const PnebiLUT& pnebiLUT, std::vector<double>& coeffs) {
        double cth2, sth2;
        //fastCosSin(2.0 * phi, cth2, sth2); 
        cth2 = cos(2.0 * phi);
        sth2 = sin(2.0 * phi);
        updateARSF2CoeffRecursDownLUT(lambda, cth2, sth2, factor, n, pnebiLUT, coeffs);
    }

    // --------------------------------------------------------
    // CLASS ISOTROPIC KERNEL 
    // --------------------------------------------------------

    class ArsKernelIsotropic2d {
    public:

        enum class ComputeMode : unsigned int {
            PNEBI_DOWNWARD = 0, PNEBI_LUT = 1
        };

        /**
         * Default constructor.
         */
        ArsKernelIsotropic2d();

        /**
         * Initializes the isotropic kernel of the given pair of identical isotropic
         * Gaussians with standard deviation sigma (i.e. covariance matrix 
         *       Cov = sigma^2 * I     where I is the identity matrix).
         * @param mean1 mean value of first Gaussian 
         * @param mean2 mean value of second Gaussian 
         * @param sigma the standard deviation of the two identical Gaussian
         */
        ArsKernelIsotropic2d(const Vector2& mean1, const Vector2& mean2, double sigma);

        /**
         * Initializes the isotropic kernel of the given pair of isotropic Gaussians
         * using their mean values and standard deviations. 
         * @param mean1 mean value of first Gaussian 
         * @param mean2 mean value of second Gaussian 
         * @param sigma1 the standard deviation of the first Gaussian
         * @param sigma2 the standard deviation of the second Gaussian
         */
        ArsKernelIsotropic2d(const Vector2& mean1, const Vector2& mean2, double sigma1, double sigma2);

        /**
         * Destructor. 
         */
        virtual ~ArsKernelIsotropic2d();

        /**
         * Initializes the isotropic kernel of the given pair of identical isotropic
         * Gaussians with standard deviation sigma (i.e. covariance matrix 
         *       Cov = sigma^2 * I     where I is the identity matrix).
         * @param mean1 mean value of first Gaussian 
         * @param mean2 mean value of second Gaussian 
         * @param sigma the standard deviation of the two identical Gaussian
         */
        void init(const Vector2& mean1, const Vector2& mean2, double sigma);

        /**
         * Initializes the isotropic kernel of the given pair of isotropic Gaussians
         * using their mean values and standard deviations. 
         * @param mean1 mean value of first Gaussian 
         * @param mean2 mean value of second Gaussian 
         * @param sigma1 the standard deviation of the first Gaussian
         * @param sigma2 the standard deviation of the second Gaussian
         */
        void init(const Vector2& mean1, const Vector2& mean2, double sigma1, double sigma2);
        
        /**
         * Initializes the LUT (look-up table) of the PNEBI (Product of Negative 
         * Exponential and Bessel I) required for fast computation of Fourier
         * coefficients of isotropic kernel. 
         * @param n the order of the PNEBI (same as the order of the truncated Fourier series)
         * @param tol the tolerance on the maximum error in PNEBI computation. 
         */
        void initPnebiLut(int n, double tol);

        /**
         * Returns the module of the sinusoidal numerator. 
         */
        double getLambdaNorm() const {
            return lambdaSqNorm_;
        }
        
        /**
         * Returns the variance associated to the ARS Kernel. 
         * @return 
         */
        double getVariance() const {
            return sigmaValSq_;
        }

        /**
         * Return the phase of the sinusoidal numerator. 
         * @return 
         */
        double getPhi() const {
            return phi_;
        }
        
        /**
         * Sets the mode for computing ARS coefficients in the case of isotropic kernels. 
         * @param mode the desired mode
         */
        void setComputeMode(ComputeMode mode) {
            mode_ = mode;
        }

        /**
         * Returns a string description of the current mode. 
         * @return 
         */
        const std::string& getComputeModeName() const {
            int idx = static_cast<int>(mode_);
            return MODE_NAME[idx];
        }

 /**
         * Computes the coefficients of Fourier series expansion of the kernel. 
         * The series is truncated to n-th order. 
         * @param nFourier maximum order of the Fourier series 
         * @param coeffs the computed coefficients 
         */
        void computeFourier(int nFourier, std::vector<double>& coeffs);
        
        /**
         * Computes the coefficients of Fourier series expansion of the kernel 
         * and ADD them to the passed ones. (It does not reset them to zero!)
         * The series is truncated to n-th order. 
         * @param nFourier maximum order of the Fourier series 
         * @param coeffs the computed coefficients 
         */
        void updateFourier(int nFourier, std::vector<double>& coeffs);
        
        /**
         * Computes the coefficients of Fourier series expansion of the kernel 
         * and ADD them to the passed ones. (It does not reset them to zero!)
         * The series is truncated to n-th order. 
         * The addeted coefficients are multiplied by w. 
         * @param nFourier maximum order of the Fourier series 
         * @param coeffs the computed coefficients 
         * @param weight the weight of this specific kernel
         */
        void updateFourier(int nFourier, std::vector<double>& coeffs, double weight);
        
        /**
         * Copmputes the value of the isotropic kernel at a given angle. 
         * @param t
         * @return 
         */
        inline double value(double t) const {
            return exp(-lambdaSqNorm_ * (1.0 + cos(2.0 * t - 2.0 * phi_))) / sqrt(2.0 * M_PI * sigmaValSq_);
        }

    private:
        static std::array<std::string, 2> const MODE_NAME;
        
        double lambdaSqNorm_;
        double sigmaValSq_;
        double phi_;
        PnebiLUT pnebiLut_;
        ComputeMode mode_;
    };

}

#endif /* ISOTROPICKERNEL_H */

