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

class NonIsotropicKernel {
public:
  /**
   * Creates a flat kernel
   */
  NonIsotropicKernel();
  
  /**
   * Constructor of non-isotropic Kernel of Angular Randon Specturm associated to two Gaussian distributions. 
   * @param mean1 mean value of first gaussian
   * @param covar1 covariance matrix of first gaussian
   * @param mean2 mean value of second gaussian
   * @param covar2 covariance matrix of second gaussian
   */
  NonIsotropicKernel(const Vector2& mean1, const Matrix2& covar1, const Vector2& mean2, const Matrix2& covar2);
  
  /**
   */
  ~NonIsotropicKernel();

  /**
   * Computes the kernel parameters and save them. 
   * @param mean1 mean value of first gaussian
   * @param covar1 covariance matrix of first gaussian
   * @param mean2 mean value of second gaussian
   * @param covar2 covariance matrix of second gaussian
   */
  void init(const Vector2& mean1, const Matrix2& covar1, const Vector2& mean2, const Matrix2& covar2);
  
  inline double value(double t) const {
      double var = sigmaMod_ * (1.0 + sigmaDif_ * cos(2.0 * t - 2.0 * sigmaAng_));
      //double mean = muMod_ * cos(2.0 * t - 2.0 * muAng_);
      return exp(-0.5 * muMod_ * muMod_ * cos(2.0 * t - 2.0 * muAng_) / var) / sqrt(2 * M_PI * var);
  }
  
  void computeFourier(int n, int k, std::vector<double>& coeffs) const;

private:
  double muMod_; 
  double muAng_;
  double sigmaMod_;
  double sigmaAng_;
  double sigmaDif_;
};

}
