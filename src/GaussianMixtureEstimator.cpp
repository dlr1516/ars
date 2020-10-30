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

    GaussianMixtureEstimator::GaussianMixtureEstimator() : means_(), covars_(), weights_() {
    }

    GaussianMixtureEstimator::~GaussianMixtureEstimator() {
    }

    //-----------------------------------------------------
    // GaussianMixtureEstimatorScan
    //-----------------------------------------------------

    GaussianMixtureEstimatorScan::GaussianMixtureEstimatorScan() : GaussianMixtureEstimator() {
    }

    GaussianMixtureEstimatorScan::~GaussianMixtureEstimatorScan() {
    }

    void GaussianMixtureEstimatorScan::compute(const VectorVector2& samples) {
        std::deque<IndexInterval> intervals;
        IndexInterval interv, interv1, interv2;
        double dist, distMax;
        int farthest;

        // Splits the scan points into intervals when a gap between consecutive 
        // points is found 
        interv.first = 0;
        for (int i = 1; i < samples.size(); ++i) {
            dist = (samples[i] - samples[i - 1]).norm();
            if (dist > distanceGap_) {
                interv.second = i - 1;
                intervals.push_back(interv);
                interv.first = i;
            }
        }
        interv.second = samples.size() - 1;
        intervals.push_back(interv);

        // Searches for aligned points in interval 
        while (!intervals.empty()) {
            // Extracts the first interval
            interv = intervals.front();
            intervals.pop_front();
            // Checks if there is a lost
            findFarthest(samples, interv.first, interv.second, farthest, distMax);
            if (distMax > distanceSplit_) {
                interv1.first = interv.first;
                interv1.second = farthest;
                interv2.first = farthest;
                interv2.second = interv.second;
            } else {

            }
        }
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
        double sigmaLowSquare = sigmaLow_ * sigmaLow_;



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
            covar << sigmaLowSquare, 0.0,
                    0.0, sigmaLowSquare;
        } else {
            covar = covar / (last - first);
            diagonalize(covar, l, v);
            if (l(0, 0) < sigmaLowSquare)
                l(0, 0) = sigmaLowSquare;
            if (l(1, 1) < sigmaLowSquare)
                l(1, 1) = sigmaLowSquare;
            covar = v * l * v.transpose();
        }
    }

}

