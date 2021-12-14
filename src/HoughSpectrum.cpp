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
#include <ars/HoughSpectrum.h>
#include <cassert>

using namespace std;

namespace ars {

    HoughSpectrum::HoughSpectrum()
    : thetaNum_(360),
    rhoNum_(1000),
    thetaStep_(0.0087266), rhoStep_(0.02),
    hough_(thetaNum_, rhoNum_),
    spectrum_(thetaNum_),
    cosLut_(thetaNum_), sinLut_(thetaNum_) {
        for (unsigned int i = 0; i < thetaNum_; ++i) {
            cosLut_(i) = cos(thetaStep_ * i);
            sinLut_(i) = sin(thetaStep_ * i);
        }
    }

    HoughSpectrum::HoughSpectrum(double thetaStep, double rhoStep, double rhoMax)
    : thetaNum_(static_cast<unsigned int> (ceil(M_PI / thetaStep))),
    rhoNum_(2 * static_cast<unsigned int> (ceil(rhoMax / rhoStep))),
    thetaStep_(thetaStep), rhoStep_(rhoStep),
    hough_(thetaNum_, rhoNum_),
    spectrum_(thetaNum_),
    cosLut_(thetaNum_), sinLut_(thetaNum_) {
        for (unsigned int i = 0; i < thetaNum_; ++i) {
            cosLut_(i) = cos(thetaStep_ * i);
            sinLut_(i) = sin(thetaStep_ * i);
        }
    }

    HoughSpectrum::~HoughSpectrum() {
    }

    void HoughSpectrum::init(double thetaStep, double rhoStep, double rhoMax) {
        thetaStep_ = thetaStep;
        rhoStep_ = rhoStep;
        thetaNum_ = (unsigned int) ceil(M_PI / thetaStep);
        rhoNum_ = (unsigned int) ceil(rhoMax / rhoStep);

        hough_.resize(thetaNum_, rhoNum_);
        spectrum_.resize(thetaNum_);
        cosLut_.resize(thetaNum_);
        sinLut_.resize(thetaNum_);

        for (unsigned int i = 0; i < thetaNum_; ++i) {
            cosLut_(i) = cos(thetaStep_ * i);
            sinLut_(i) = sin(thetaStep_ * i);
        }
    }

    const Eigen::MatrixXd& HoughSpectrum::hough() const {
        return hough_;
    }

    double HoughSpectrum::hough(double theta, double rho) const {
        int ith = thetaToIdx(theta);
        int irh = rhoToIdx(rho);
        if (0 <= ith && ith < hough_.rows() && 0 <= irh && irh < hough_.cols()) {
            return hough_(ith, irh);
        }
        return 0.0;
    }

    const Eigen::VectorXd& HoughSpectrum::spectrum() const {
        return spectrum_;
    }

    const double HoughSpectrum::spectrum(double theta) const {
        int ith = thetaToIdx(theta);
        if (0 <= ith && ith < hough_.rows()) {
            return spectrum_(ith);
        }
        return 0.0;
    }

    double HoughSpectrum::idxToRho(int idx) const {
        return (rhoStep_ * (idx - rhoNum_));
    }

    int HoughSpectrum::rhoToIdx(double rho) const {
        return ((int) round(rho / rhoStep_) + rhoNum_ / 2);
    }

    int HoughSpectrum::thetaToIdx(double theta) const {
        int idx = (int) round(theta / thetaStep_);
        int thetaNum2 = 2 * thetaNum_;
        idx = ((idx % thetaNum2) + thetaNum2) % thetaNum2;
        return idx;
    }

} // end of namespace

