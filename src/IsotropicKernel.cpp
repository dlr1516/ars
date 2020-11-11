/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   IsotropicKernel.cpp
 * Author: pmicelli
 * 
 * Created on November 9, 2020, 11:59 AM
 */

#include <ars/IsotropicKernel.h>

namespace ars {

    IsotropicKernel::IsotropicKernel() : lambdaSqNorm_(0.0), sigmaValSq_(1.0), phi_(0.0), pnebiLut_(40, 0.001), mode_(ComputeMode::PNEBI_LUT) {
    }

    IsotropicKernel::IsotropicKernel(const Vector2& mean1, const Vector2& mean2, double sigma)
    : pnebiLut_(40, 0.001), mode_(ComputeMode::PNEBI_LUT) {
        init(mean1, mean2, sigma);
    }

    IsotropicKernel::IsotropicKernel(const Vector2& mean1, const Vector2& mean2, double sigma1, double sigma2)
    : pnebiLut_(40, 0.001), mode_(ComputeMode::PNEBI_LUT) {
        init(mean1, mean2, sigma1, sigma2);
    }

    IsotropicKernel::~IsotropicKernel() {
    }

    void IsotropicKernel::init(const Vector2& mean1, const Vector2& mean2, double sigma) {
        init(mean1, mean2, sigma, sigma);
    }

    void IsotropicKernel::init(const Vector2& mean1, const Vector2& mean2, double sigma1, double sigma2) {
        double dx, dy;

        dx = mean2(0) - mean1(0);
        dy = mean2(1) - mean1(1);
        phi_ = atan2(dy, dx);
        sigmaValSq_ = sigma1 * sigma1 + sigma2 * sigma2;
        lambdaSqNorm_ = 0.25 * (dx * dx + dy * dy) / sigmaValSq_;
    }

    void IsotropicKernel::initPnebiLut(int n, double tol) {
        pnebiLut_.init(n, tol);
    }

    void IsotropicKernel::computeFourier(int nFourier, std::vector<double>& coeffs) const {
        if (coeffs.size() != 2 * nFourier + 2) {
            coeffs.resize(2 * nFourier + 2);
        }
        std::fill(coeffs.begin(), coeffs.end(), 0.0);

        updateFourier(nFourier, coeffs);
    }

    void IsotropicKernel::updateFourier(int nFourier, std::vector<double>& coeffs) const {
        double w = 1.0 / sqrt(2.0 * M_PI * sigmaValSq_);
        
        if (coeffs.size() != 2 * nFourier + 2) {
            coeffs.resize(2 * nFourier + 2);
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

    std::array<std::string, 2> const IsotropicKernel::MODE_NAME{"PNEBI_DOWNWARD", "PNEBI_LUT"};


} // end of namespace

