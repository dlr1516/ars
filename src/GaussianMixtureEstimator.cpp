/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   GaussianMixtureEstimator.cpp
 * Author: pmicelli
 * 
 * Created on October 30, 2020, 8:33 AM
 */

#include <ars/GaussianMixtureEstimator.h>
#include <vector>
#include <deque>

#include "ars/utils.h"

namespace ars {

    //-----------------------------------------------------
    // GaussianMixtureEstimator
    //-----------------------------------------------------

    GaussianMixtureEstimator::GaussianMixtureEstimator() : gaussians_() { //: means_(), covars_(), weights_() {
    }

    GaussianMixtureEstimator::~GaussianMixtureEstimator() {
    }

    //-----------------------------------------------------
    // GaussianMixtureEstimatorScan
    //-----------------------------------------------------

    GaussianMixtureEstimatorScan::GaussianMixtureEstimatorScan() 
     : GaussianMixtureEstimator(), intervals_(), distanceGap_(0.8), distanceSplit_(0.5), sigmaMin_(0.1) {
    }

    GaussianMixtureEstimatorScan::~GaussianMixtureEstimatorScan() {
    }

    void GaussianMixtureEstimatorScan::compute(const VectorVector2& samples) {
        Vector2 mean;
        Matrix2 covar;
        std::deque<IndexInterval> intervals;
        IndexInterval interv, interv1, interv2;
        double dist, distMax, w;
        int farthest;
        
        int sum = 0;

        // Splits the scan points into intervals when a gap between consecutive 
        // points is found 
        interv.first = 0;
        for (int i = 1; i < samples.size(); ++i) {
            dist = (samples[i] - samples[i - 1]).norm();
            if (dist > distanceGap_) {
                interv.last = i - 1;
                interv.num = interv.last - interv.first + 1; 
                intervals.push_back(interv);
                ARS_PRINT("interv [" << interv.first << ", " << interv.last << "] num " << interv.num);
                interv.first = i;
            }
        }
        interv.last = samples.size() - 1;
        interv.num = interv.last - interv.first + 1; 
        intervals.push_back(interv);
        
        std::cout << "\n----\n" << std::endl;

        // Searches for aligned points in interval 
        while (!intervals.empty()) {
            // Extracts the first interval
            interv = intervals.front();
            intervals.pop_front();
            // Checks if the interval is split according to a policy based on 
            // distance of farthest point from segment
            findFarthest(samples, interv.first, interv.last, farthest, distMax);
            if (distMax > distanceSplit_) {
                // Interval is split at farthest point. Formally:
                // - interv1: [interv.first, farthest-1] (farthest NOT included**)
                // - interv2: [farthest, interv.last] 
                // ** the fathest is not included in the first interval, but it's used 
                //    in the computation of the gaussian!!! So we have an overlap:
                // - interv1: [interv.first, farthest] (farthest NOT included**)
                interv1.first = interv.first;
                interv1.last = farthest;        
                interv1.num = farthest - interv.first;
                interv2.first = farthest;
                interv2.last = interv.last;
                interv2.num = interv.num - interv1.num;
                intervals.push_front(interv1);
                intervals.push_front(interv2);
            } else {
                ARS_PRINT("interv [" << interv.first << ", " << interv.last << "] num " << interv.num);
                Gaussian g;
                estimateGaussian(samples, interv.first, interv.last, g.mean, g.covar);
                w = interv.num * 1.0 / samples.size();                
                g.weight = w;
                gaussians_.push_front(g);
                //means_.push_back(mean);
                //covars_.push_back(covar);
                //weights_.push_back(w);
                intervals_.push_front(interv);
                sum += interv.num;
            }
        }
        
        ARS_VARIABLE(sum);
    }

    void GaussianMixtureEstimatorScan::findFarthest(const VectorVector2& points, int first, int last, int& farthest, double& distMax) const {
        Vector2 dp;
        double t, ct, st, r, dist;

        ARS_ASSERT(0 <= first && last < points.size() && first <= last);
        dp = points[last] - points[first];
        t = atan2(dp(0), -dp(1));
        ct = cos(t);
        st = sin(t);
        r = points[first](0) * ct + points[first](1) * st;

        distMax = 0.0;
        for (int i = first; i <= last; ++i) {
            dist = fabs(points[i](0) * ct + points[i](1) * st - r);
            if (dist > distMax) {
                distMax = dist;
                farthest = i;
            }
        }
    }

    void GaussianMixtureEstimatorScan::estimateGaussian(const VectorVector2& points, int first, int last, Vector2& mean, Matrix2& covar) const {
        Matrix2 l, v;
        Vector2 tmp;
        double sigmaMinSquare = sigmaMin_ * sigmaMin_;

        ARS_ASSERT(first >= 0 && last < points.size());

        // Computes the mean value vector
        mean = Vector2::Zero();
        for (int i = first; i <= last; ++i) {
            mean += points[i];
        }
        mean = mean / (last - first + 1);

        // Computes the covariance
        covar = Matrix2::Zero();
        for (int i = first; i <= last; ++i) {
            tmp = (points[i] - mean);
            covar = tmp * tmp.transpose();
        }

        if (first == last) {
            // Only one point: use the point uncertainty
            covar << sigmaMinSquare, 0.0,
                    0.0, sigmaMinSquare;
        } else {
            covar = covar / (last - first);
            diagonalize(covar, l, v);
            if (l(0, 0) < sigmaMinSquare)
                l(0, 0) = sigmaMinSquare;
            if (l(1, 1) < sigmaMinSquare)
                l(1, 1) = sigmaMinSquare;
            covar = v * l * v.transpose();
        }
    }

}

