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
#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <ars/definitions.h>
#include <ars/functions.h>
#include <ars/BBOptimizer1d.h>

namespace ars {

    // --------------------------------------------------------
    // ARS-2D FUNCTIONS
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

    /** Computes coeffients of Fourier series Correlation. 
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

        enum ComputeMode {
            PNEBI_DOWNWARD, PNEBI_LUT
        };

        /** Default constructor. 
         */
        AngularRadonSpectrum2d();

        /** Destructor. 
         */
        ~AngularRadonSpectrum2d();

        /** Sets the order of Fourier approximation of ARS.
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

        void setComputeMode(ComputeMode mode) {
            mode_ = mode;
        }

        const std::string& getComputeModeName() const {
            return MODE_NAME[mode_];
        }

        /** Returns const reference to ARSF coefficients obtained from downward recursion. 
         */
        const std::vector<double>& coefficients() const {
            return coeffs_;
        }

        /** Inserts the given points and computes all the data about point pairs.
         */
        void insertIsotropicGaussians(const Vector2Vector& means, double sigma);

        /** Inserts the given points and computes all the data about point pairs.
         */
        void insertIsotropicGaussians(const Vector2Vector& means, const std::vector<double>& sigmas);

        /** Initializes LUT (the LUT is used by initARSFRecursDownLUT).
         */
        void initLUT(double precision = 0.001) {
            pnebiLut_.init(arsfOrder_, precision);
        }

        /** Evaluates the ARSF using the coefficients obtained from downward recursion. 
         */
        double eval(double theta) const;

        /** Finds the maximum of ARSF.
         */
        double findMax() const {
            double arsfMax, thetaMax;
            findGlobalMaxBBFourier(coeffs_, 0, M_PI, thetaToll_, 10.0, thetaMax, arsfMax);
            return arsfMax;
        }

        /** Finds the maximum of ARSF.
         */
        double findMax(double& thetaMax) const {
            double arsfMax;
            findGlobalMaxBBFourier(coeffs_, 0, M_PI, thetaToll_, 10.0, thetaMax, arsfMax);
            return arsfMax;
        }

    protected:
        static std::array<std::string, 2> const MODE_NAME;

        std::vector<double> coeffs_;
        int arsfOrder_;
        double thetaToll_;
        int threadNumOMP_;
        // PNEBI LUT
        PnebiLUT pnebiLut_;
        ComputeMode mode_;
    };

} // end of namespace

