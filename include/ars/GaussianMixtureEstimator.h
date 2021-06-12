/**
 * ARS - Angular Radon Spectrum 
 * Copyright (C) 2017 Dario Lodi Rizzini.
 *           (C) 2021 Dario Lodi Rizzini, Ernesto Fontana.
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
#include <deque>
#include <Eigen/Dense>
#include <ars/definitions.h>
//#include <ars/MortonOrderedPoints.h>
#include <ars/DisjointSet.h>

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

        struct Gaussian {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            Vector2 mean;
            Matrix2 covar;
            double weight;

            double eval(const Vector2& v) const {
                double k = 1.0 / sqrt(2.0 * M_PI * covar.determinant());
                double arg = (v - mean).transpose() * covar.inverse() * (v - mean);
                return k * exp(-0.5 * arg);
            }
        };
        using VectorGaussian = std::deque<Gaussian, Eigen::aligned_allocator<Gaussian> >;

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
            //return weights_.size();
            return gaussians_.size();
        }

        /**
         * Returns the estimated mean value of i-th Gaussian distribution in the mixture.
         * @param i the index of the distribution/hypothesis
         * @return the mean vector
         */
        const Vector2& mean(int i) const {
            //            ARS_ASSERT(0 <= i && i < means_.size());
            //            return means_.at(i);
            ARS_ASSERT(0 <= i && i < gaussians_.size());
            return gaussians_[i].mean;
        }

        /**
         * Returns the estimated covariance of i-th Gaussian distribution in the mixture.
         * @param i the index of the distribution/hypothesis
         * @return the covariance matrix
         */
        const Matrix2& covariance(int i) const {
            //            ARS_ASSERT(0 <= i && i < covars_.size());
            //            return covars_.at(i);
            ARS_ASSERT(0 <= i && i < gaussians_.size());
            return gaussians_[i].covar;
        }

        /**
         * Returns the estimated weight of i-th Gaussian distribution in the mixture,
         * i.e. the probability that the i-th component/hypothesis is drawn. 
         * @param i the index of the distribution/hypothesis
         * @return the weight
         */
        double weight(int i) const {
            //            ARS_ASSERT(0 <= i && i < weights_.size());
            //            return weights_.at(i);
            ARS_ASSERT(0 <= i && i < gaussians_.size());
            return gaussians_[i].weight;
        }

        /**
         * Exports the Gaussian mixture parameters, i.e. means, covariances and weights, 
         * into separate vectors. 
         * @param means std::vector of mean vectors
         * @param covariances std::vector of covariance matrices
         * @param weights std::vector of weights
         */
        void exportGaussians(VectorVector2& means, VectorMatrix2& covariances, std::vector<double>& weights) const;

        /**
         * Executes Expectation Maximization (EM) updating the Gaussian 
         * mixture stored in the class.  
         * @param samples vector of samples
         * @param stepNum number of iteration of EM
         */
        void executeEM(const VectorVector2& samples, int stepNum = 1);

    protected:
        //        VectorVector2 means_;
        //        VectorMatrix2 covars_;
        //        std::vector<double> weights_;
        VectorGaussian gaussians_;
    };

    //-----------------------------------------------------
    // GaussianMixtureEstimatorScan
    //-----------------------------------------------------

    class GaussianMixtureEstimatorScan : public GaussianMixtureEstimator {
    public:

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        //using IndexInterval = std::pair<int, int>;
        struct IndexInterval {
            int first;
            int last;
            int num;
        };

        /**
         * Default constructor. 
         */
        GaussianMixtureEstimatorScan();

        /**
         * Destructor. 
         */
        virtual ~GaussianMixtureEstimatorScan();

        /**
         * Sets the threshold above which a gap between consecutive points is detected. 
         * @param dg the distance gap threshold
         */
        void setDistanceGap(double dg) {
            distanceGap_ = dg;
        }

        /**
         * Sets the splitting distance for segment detection
         * @param ds the distance threshold to split 
         */
        void setDistanceSplit(double ds) {
            distanceSplit_ = ds;
        }

        /**
         * Sets the minimum value of standard deviation of Gaussians.
         * @param sm the minimum standard deviation 
         */
        void setSigmaMin(double sm) {
            sigmaMin_ = sm;
        }

        /**
         * Computes the Gaussian parameters from the given samples.
         * @param samples sorted in counter-clockwise order
         */
        virtual void compute(const VectorVector2& samples);

        /**
         * Returns the i-th interval. 
         * @param i
         */
        const IndexInterval& interval(int i) const {
            ARS_ASSERT(0 <= i && i < intervals_.size());
            return intervals_.at(i);
        }

    private:
        std::deque<IndexInterval> intervals_; // used for debug 
        double distanceGap_;
        double distanceSplit_;
        double sigmaMin_;

        /**
         * Finds the farthest point in the set from the line through points first and last,
         * i.e. points[first] and points[last]. 
         * The farthest index belongs to interval first and last. 
         * @param points the complete set of points
         * @param first the index of the first point in the segment interval
         * @param last the index of the last point in the segment interval
         * @param farthest the index of the farthest point from the line
         * @param distMax the distance of the farthest point from the line
         */
        void findFarthest(const VectorVector2& points, int first, int last, int& farthest, double& distMax) const;

        /**
         * Computes the Gaussian mean and covariance matrix of points in interval
         * between first and last. 
         * @param points the set of points
         * @param first the index of the first point 
         * @param last the intex of the last point
         * @param mean the mean vector
         * @param covar the covariance matrix
         */
        void estimateGaussianFromPoints(const VectorVector2& points, int first, int last, Vector2& mean, Matrix2& covar) const;

        /**
         * Computes the Gaussian distribution, i.e. its parameters, assuming the input points 
         * are not samples of the Gaussian, but rather approximately uniform sampled points
         * on the segment. 
         * @param points the set of points
         * @param first the index of the first point 
         * @param last the intex of the last point
         * @param mean the mean vector
         * @param covar the covariance matrix
         */
        void estimateGaussianFromSegment(const VectorVector2& points, int first, int last, Vector2& mean, Matrix2& covar) const;

    };

    //-----------------------------------------------------
    // GaussianMixtureEstimatorMeanShift
    //-----------------------------------------------------

    class GaussianMixtureEstimatorMeanShift : public GaussianMixtureEstimator {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        /**
         * Default constructor.
         */
        GaussianMixtureEstimatorMeanShift();

        /**
         * Destructor. 
         */
        virtual ~GaussianMixtureEstimatorMeanShift();

        //            /**
        //             * Sets the number of Gaussian kernels of the mixture. 
        //             * @param knum
        //             */
        //            void setKernelNum(int knum) {
        //                kernelNum_ = knum;
        //            }

        /**
         * Sets the minimum value of standard deviation of Gaussians.
         * @param sm the minimum standard deviation 
         */
        void setSigmaMin(double sm) {
            sigmaMin_ = sm;
        }

        /**
         * Sets the cluster threshold to create a new cluster. The distances are 
         * normalized. For example, given two points p and q, their distance is:
         * 
         *   dist(p, q) = |p - q| / sigmaMin
         * 
         * @param cd the value of cluster distance
         */
        void setClusterDistance(double cd) {
            clusterDist_ = cd;
        }
        
        /**
         * Sets the tolerance on the distance between the shifted means inside 
         * the same cluster. 
         * The clustering algorithm is stopped when all the intra-cluster distances
         * are less than meanShiftTol_. 
         * @param mst
         */
        void setMeanShiftTol(double mst) {
            meanShiftTol_ = mst;
        }
        
        /**
         * Sets the maximum number of iterations of mean-shift clustering algorithm. 
         * @param inmax
         */
        void setIterationNumMax(int inmax) {
            iterationNumMax_ = inmax;
        }

        /**
         * Computes the Gaussian parameters from the given samples.
         * @param samples sorted in counter-clockwise order
         */
        virtual void compute(const VectorVector2& samples);


    private:
        int kernelNum_;
        double sigmaMin_;
        double clusterDist_;
        double meanShiftTol_;
        int iterationNumMax_;

        void updateMeans(const VectorVector2& meansCurr, VectorVector2& meansNext, DisjointSet& clusterLabels, std::vector<double>& clusterIntraDistMax) const;
    };

    //-----------------------------------------------------
    // GaussianMixtureEstimatorHierarchical
    //-----------------------------------------------------

    /**
     * Computes a Gaussian Mixture approximately exploiting some of the ideas in 
     * 
     * Benjamin Eckart, Kihwan Kim Jan Kau,
     * "HGMR: Hierarchical Gaussian Mixtures for Adaptive 3D Registration",
     * ECCV 2018. 
     * 
     */
    //    class GaussianMixtureEstimatorHierarchical : public GaussianMixtureEstimator {
    //    public:
    //        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    //
    //                // Private
    //                using PointContainer = MortonOrderedPoints2d16h;
    //        using Iterator = PointContainer::Iterator;
    //        using ConstIterator = PointContainer::ConstIterator;
    //
    //        /**
    //         * Default constructor.
    //         */
    //        GaussianMixtureEstimatorHierarchical();
    //
    //        /**
    //         * Destructor. 
    //         */
    //        virtual ~GaussianMixtureEstimatorHierarchical();
    //
    //        /**
    //         * Sets the minimum value of standard deviation of Gaussians.
    //         * @param sm the minimum standard deviation 
    //         */
    //        void setSigmaMin(double sm) {
    //            sigmaMin_ = sm;
    //        }
    //
    //        /**
    //         * Computes the Gaussian parameters from the given samples.
    //         * @param samples sorted in counter-clockwise order
    //         */
    //        virtual void compute(const VectorVector2& samples);
    //
    //    private:
    //        PointContainer data_;
    //        double sigmaMin_;
    //
    //        void estimateGaussianFromPoints(const ConstIterator& beg, const ConstIterator& end, Vector2& mean, Matrix2& covar) const;
    //    };

} // end of namespace

#endif /* GAUSSIANMIXTUREESTIMATOR_H */

