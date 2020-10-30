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
#ifndef GAUSSIANMIXTUREESTIMATOR_H
#define GAUSSIANMIXTUREESTIMATOR_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <ars/definitions.h>

namespace ars {
    
    //-----------------------------------------------------
    // GaussianMixtureEstimator
    //-----------------------------------------------------

    /**
     * Class GaussianMixtureEstimator provides a general interface for the estimators 
     * of Gaussian Mixture Models (GMMs) that compute the Gaussian parameters 
     * (mean vectors, covariance matrices, weights, etc.) from observed samples. 
     */
    class GaussianMixtureEstimator {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        /**
         * Default constructor.
         */
        GaussianMixtureEstimator();
        
        /**
         * Destructor. 
         */
        virtual ~GaussianMixtureEstimator();
        
        /**
         * Computes the Gaussian parameters from the given samples.
         * @param samples
         */
        virtual void compute(const VectorVector2& samples) = 0;
        
        /**
         * Returns the number of components/hypotheses of the mixture.
         * @return 
         */
        int size() const {
            return weights_.size();
        }
        
        /**
         * Returns the estimated mean value of i-th Gaussian distribution in the mixture.
         * @param i the index of the distribution/hypothesis
         * @return the mean vector
         */
        const Vector2& mean(int i) const {
            ARS_ASSERT(0 <= i && i < means_.size());
            return means_.at(i);
        }
        
        /**
         * Returns the estimated covariance of i-th Gaussian distribution in the mixture.
         * @param i the index of the distribution/hypothesis
         * @return the covariance matrix
         */
        const Matrix2& covariance(int i) const {
            ARS_ASSERT(0 <= i && i < covars_.size());
            return covars_.at(i);
        }
        
        /**
         * Returns the estimated weight of i-th Gaussian distribution in the mixture,
         * i.e. the probability that the i-th component/hypothesis is drawn. 
         * @param i the index of the distribution/hypothesis
         * @return the weight
         */
        double weight(int i) const {
            ARS_ASSERT(0 <= i && i < weights_.size());
            return weights_.at(i);
        }
        
    protected:
        VectorVector2 means_;
        VectorMatrix2 covars_;
        std::vector<double> weights_;
    };
    
    //-----------------------------------------------------
    // GaussianMixtureEstimatorScan
    //-----------------------------------------------------
    
    class GaussianMixtureEstimatorScan : public GaussianMixtureEstimator {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        using IndexInterval = std::pair<int, int>;
        
        /**
         * Default constructor. 
         */
        GaussianMixtureEstimatorScan();
        
        /**
         * Destructor. 
         */
        virtual ~GaussianMixtureEstimatorScan();
        
        /**
         * Computes the Gaussian parameters from the given samples.
         * @param samples sorted in counter-clockwise order
         */
        virtual void compute(const VectorVector2& samples);
        
    private:
        double distanceGap_;
        double distanceSplit_;
        double sigmaLow_;
        
        void findFarthest(const VectorVector2& points, int first, int last, int& farthest, double& distMax) const;
        
        void estimateGaussian(const VectorVector2& points, int first, int last, Vector2& mean, Matrix2& covar) const;
    };

} // end of namespace

#endif /* GAUSSIANMIXTUREESTIMATOR_H */

