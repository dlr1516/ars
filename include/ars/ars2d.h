/**
 * ARS - Angular Radon Spectrum 
 * Copyright (C) 2017-2020 Dario Lodi Rizzini.
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
#include <ars/functions.h>
#include <ars/BBOptimizer1d.h>
#include <ars/ArsKernelIsotropic2d.h>
#include <ars/ArsKernelAnisotropic2d.h>

namespace ars {

    /** Computes coefficients of Fourier series Correlation. 
     * Given two series Fsrc(t) and Fdst(t) with coefficients:
     *   Fsrc(t) = \sum_{i=0}^{n} ( fourierSrc[2*i] * cos(2*i*t) + fourierSrc[2*i+1] * sin(2*i*t) )
     *   Fdst(t) = \sum_{i=0}^{n} ( fourierDst[2*i] * cos(2*i*t) + fourierDst[2*i+1] * sin(2*i*t) )
     * the correlation of the two function is the integral over period T:
     *   Fcor(h) = \int_{0}^{T} Fsrc(t+h) Fdst(t) dt
     * where Fcor(h) is still represented as a Fourier serie with coeffients fourierCor[]. 
     */
    void computeFourierCorr(const std::vector<double>& fourierSrc, const std::vector<double>& fourierDst, std::vector<double>& fourierCor);

    // --------------------------------------------------------
    // ARS 2D CLASS
    // --------------------------------------------------------

    /** Class for computing Angular Radon Spectum (ARS) of a set of points.
     */
    class AngularRadonSpectrum2d {
    public:

//        enum ComputeMode {
//            PNEBI_DOWNWARD, PNEBI_LUT
//        };

        /** Default constructor. 
         */
        AngularRadonSpectrum2d();
        
        /** Default constructor. 
         */
        AngularRadonSpectrum2d(const std::vector<double>& coeffs);

        /** Destructor. 
         */
        ~AngularRadonSpectrum2d();

        /** Sets the order of truncated Fourier series expansion of ARS.
         * WARNING: It resets the coefficients!
         */
        void setARSFOrder(int n) {
            arsfOrder_ = n;
            coeffs_.resize(2 * n + 2);
            //initLUT();
        }

        /** Sets the maximum tollerance on theta during the computation of maximum.
         */
        void setThetaToll(double thetaToll) {
            thetaToll_ = thetaToll;
        }

        /** Sets the number of thread used by OpenMP routines. 
         */
        void setThreadNumOMP(int tno) {
            threadNumOMP_ = tno;
        }

        /**
         * Sets the mode for computing ARS coefficients in the case of isotropic kernels. 
         * @param mode the desired mode
         */
        void setComputeMode(ArsKernelIsotropic2d::ComputeMode mode) {
            //mode_ = mode;
            isotropicKer_.setComputeMode(mode);
        }

        /**
         * Returns a string description of the current mode. 
         * @return 
         */
        const std::string& getComputeModeName() const {
            //return MODE_NAME[mode_];
            return isotropicKer_.getComputeModeName();
        }
        
        /**
         * Sets the number of intervals used in the computation of Fourier coeffcients 
         * of anisotropic kernels. 
         * Since the closed-form equation of Fourier coefficients is unknown in 
         * anisotropic case, the numerical integration with step M_PI / anisotropicStep_
         * is used instead. 
         * @param as
         */
        void setAnisotropicStep(int as) {
            anisotropicStep_ = as;
        }

        /** Returns const reference to ARS Fourier coefficients. 
         * Coefficients are obtained from a Gaussian Mixture Model (GMM) representing
         * a point set distribution with uncertainty. 
         */
        const std::vector<double>& coefficients() const {
            return coeffs_;
        }
        
        /**
         * Sets the ARS Fourier coefficients. 
         * @warning This method is a "backdoor" w.r.t. the insert methods that 
         * computes the coefficients directly from the point set in GMM form. 
         * Accordingly it should be used carefully to avoid inconsistencies. 
         * @param coeffs the coefficients
         */
        void setCoefficients(const std::vector<double>& coeffs) {
            coeffs_ = coeffs;
        }

        /** * Inserts the given points and computes all the data about point pairs and 
         * computes the coefficients of the Fourier series representing the ARS 
         * of the point set.
         * All the Gaussian distributions are isotropic, have the same standard deviation 
         * and the same weight in the mixture. 
         * @param means mean values of the distributions
         * @param sigma the standard deviation (not variance!) of the identical isotropic distributions
         */
        void insertIsotropicGaussians(const VectorVector2& means, double sigma);

        /**
         * Inserts the given points and computes all the data about point pairs and 
         * computes the coefficients of the Fourier series representing the ARS 
         * of the point set.
         * Hypothesis: the weights of the input Gaussian distributions of the mixture
         * are assumed to be equal. 
         * @param means mean values of the distributions
         * @param sigmas standard deviations (not variances!) of the isotropic distributions
         */
        void insertIsotropicGaussians(const VectorVector2& means, const std::vector<double>& sigmas);
        
        /**
         * Inserts the given points and computes all the data about point pairs and 
         * computes the coefficients of the Fourier series representing the ARS 
         * of the point set.
         * @param means mean values of the distributions
         * @param sigmas standard deviations (not variances!) of the isotropic distributions
         * @param weights the weights of each distribution of the mixture
         */
        void insertIsotropicGaussians(const VectorVector2& means, const std::vector<double>& sigmas, const std::vector<double>& weights);
        
        /**
         * Inserts the given anisotropic gaussians. 
         * Pre-condition: means.size() == covars.size(). 
         * @param means the mean values of Gaussians PDF representing points
         * @param covars the covariance matrices of Gaussians PDF representing point uncertainties
         */
        void insertAnisotropicGaussians(const VectorVector2& means, const VectorMatrix2& covars, const std::vector<double>& weights);

        /** Initializes LUT (the LUT is used by initARSFRecursDownLUT).
         */
        void initLUT(double precision = 0.001) {
            //pnebiLut_.init(arsfOrder_, precision);
            isotropicKer_.initPnebiLut(arsfOrder_, precision);
        }

        /** Evaluates the ARS Fourier using the coefficients obtained from downward recursion. 
         */
        double eval(double theta) const;

        /** Finds the maximum of ARS Fourier.
         */
        double findMax() const {
            double arsfMax, thetaMax;
            findGlobalMaxBBFourier(coeffs_, 0, M_PI, thetaToll_, 10.0, thetaMax, arsfMax);
            return arsfMax;
        }

        /** Finds the maximum of ARS Fourier.
         */
        double findMax(double& thetaMax) const {
            double arsfMax;
            findGlobalMaxBBFourier(coeffs_, 0, M_PI, thetaToll_, 10.0, thetaMax, arsfMax);
            return arsfMax;
        }
        
        /**
         * Finds the maximum of ARS Fourier on the given interval [thetaLow, thetaUpp].
         * @param thetaLow
         * @param thetaUpp
         * @return the maximum value of correlation function. 
         */
        double findMax(double thetaLow, double thetaUpp) const {
            double arsfMax, thetaMax;
            findGlobalMaxBBFourier(coeffs_, thetaLow, thetaUpp, thetaToll_, 10.0, thetaMax, arsfMax);
            return arsfMax;
        }
        
        /**
         *  Finds the maximum of ARS Fourier on the given interval.
         * @param thetaOpt
         * @param thetaMin
         * @param thetaMax
         * @return 
         */
        double findMax(double& thetaMax, double thetaLow, double thetaUpp) const {
            double arsfMax;
            findGlobalMaxBBFourier(coeffs_, thetaLow, thetaUpp, thetaToll_, 10.0, thetaMax, arsfMax);
            return arsfMax;
        }

    protected:
        std::vector<double> coeffs_;
        ArsKernelIsotropic2d isotropicKer_;
        ArsKernelAnisotropic2d anisotropicKer_;
        int arsfOrder_;
        double thetaToll_;
        int threadNumOMP_;
        // Parameters for computation of the Fourier coefficients of anisotropic kernels
        int anisotropicStep_;
    };

} // end of namespace

