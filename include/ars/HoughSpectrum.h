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

namespace ars {

    /** Hough transform and spectrum. Lines are represented with Hessian (aka 
     * polar) parameters:
     * 
     *   x * cos(theta) + y * sin(theta) = rho
     * 
     * See the following references:
     * 
     * - A. Censi, L. Iocchi, G. Grisetti, "Scan Matching in the Hough Domain",
     *   p. 2739-2744, Proc. of IEEE ICRA 2015. 
     * - D. Lodi Rizzini, "Angular Radon Spectrum for Rotation Estimation",
     *   Pattern Recognition, Vol. 84, Dec. 2018, Pages 182-196, 
     *   DOI 10.1016/j.patcog.2018.07.017.
     * 
     */
    class HoughSpectrum {
    public:

        /** Constructor with deafult parameters. 
         */
        HoughSpectrum();

        /**
         * Constructor with initialization parameters
         * @param thetaStep the dimension of angle bin
         * @param rhoStep the dimension of range bin
         * @param rhoMax the maximum value of polar range of lines
         */
        HoughSpectrum(double thetaStep, double rhoStep, double rhoMax);

        /** Default destructor. 
         */
        virtual ~HoughSpectrum();

        /**
         * Inits params of histogram representing the discretized domain of 
         * Hough transform and spectrum. 
         * @param thetaStep the dimension of angle bin
         * @param rhoStep the dimension of range bin
         * @param rhoMax the maximum value of polar range of lines
         */
        void init(double thetaStep, double rhoStep, double rhoMax);

        /** Inserts the points and computes Hough Transform and Spectrum. 
         */
        template <typename It>
        void insertPoint(It pbeg, It pend);

        /** Returns the Hough Transform. 
         */
        const Eigen::MatrixXd& hough() const;

        /** Returns the value of Hough Transform for a specific value of theta and rho.
         * If the theta and rho are not in the domain, then it return 0.0.
         */
        double hough(double theta, double rho) const;

        /** Returns the spectrum.
         */
        const Eigen::VectorXd& spectrum() const;

        /** Returns the value of spectrum at given angle.
         * @param theta the angle value
         */
        const double spectrum(double theta) const;

    private:
        int thetaNum_;
        int rhoNum_;
        double thetaStep_;
        double rhoStep_;
        // Hough transform and spectra should be integer types. 
        // However, since their value may be very large, double type is used instead. 
        Eigen::MatrixXd hough_;
        Eigen::VectorXd spectrum_;
        Eigen::VectorXd cosLut_;
        Eigen::VectorXd sinLut_;

        double idxToRho(int idx) const;

        int rhoToIdx(double rho) const;

        int thetaToIdx(double theta) const;
    };

    // -------------------------------------------------------------
    // TEMPLATE METHOD
    // -------------------------------------------------------------

    template <typename It>
    void HoughSpectrum::insertPoint(It pbeg, It pend) {
        double rho;
        int irho;

        hough_.fill(0.0);
        spectrum_.fill(0.0);
        // for (int r = 0; r < hough_.rows(); ++r) {
        //     for (int c = 0; c < hough_.cols(); ++c) {
        //         hough_(r, c) = 0.0;
        //     }
        // }
        // for (int r = 0; r < spectrum_.rows(); ++r) {
        //     spectrum_(r) = 0.0;
        // }
        // Computes Hough
        for (It pit = pbeg; pit != pend; ++pit) {
            for (int i = 0; i < thetaNum_; ++i) {
                rho = pit->x() * cosLut_(i) + pit->y() * sinLut_(i);
                irho = rhoToIdx(rho);
                if (0 <= irho && irho < rhoNum_) {
                    hough_(i, irho) = hough_(i, irho) + 1;
                }
                //      else {
                //        std::cerr << "Out-of-bound: rho " << rho << ", theta " << (180.0 / M_PI * thetaStep_ * i) << ", "
                //          << "point [" << pit->x() << " " << pit->y() << "]\n";
                //      }
            }
        }
        // Computes theta spectrum
        //std::cout << __FILE__ << "," << __LINE__ << std::endl;
        spectrum_ = (hough_.array() * hough_.array()).rowwise().sum();
    }

} // end of namespace emotion

