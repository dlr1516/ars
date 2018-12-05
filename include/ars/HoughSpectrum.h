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
#include <Eigen/Dense>

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
        const Eigen::MatrixXd& hough() const {
            return hough_;
        }

        /** Returns the value of Hough Transform for a specific value of theta and rho.
         * If the theta and rho are not in the domain, then it return 0.0.
         */
        double hough(double theta, double rho) const {
            int ith = thetaToIdx(theta);
            int irh = rhoToIdx(rho);
            if (0 <= ith && ith < hough_.rows() && 0 <= irh && irh < hough_.cols()) {
                return hough_(ith, irh);
            }
            return 0.0;
        }

        /** Returns the spectrum.
         */
        const Eigen::VectorXd& spectrum() const {
            return spectrum_;
        }

        /** Returns the value of spectrum at given angle.
         * @param theta the angle value
         */
        const double spectrum(double theta) const {
            int ith = thetaToIdx(theta);
            if (0 <= ith && ith < hough_.rows()) {
                return spectrum_(ith);
            }
            return 0.0;
        }

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

        double idxToRho(int idx) const {
            return (rhoStep_ * (idx - rhoNum_));
        }

        int rhoToIdx(double rho) const {
            return ((int) round(rho / rhoStep_) + rhoNum_ / 2);
        }

        int thetaToIdx(double theta) const {
            int idx = (int) round(theta / thetaStep_);
            int thetaNum2 = 2 * thetaNum_;
            idx = ((idx % thetaNum2) + thetaNum2) % thetaNum2;
            return idx;
        }
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

