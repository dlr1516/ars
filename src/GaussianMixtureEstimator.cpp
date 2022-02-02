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
#include <ars/GaussianMixtureEstimator.h>
#include <vector>
#include <deque>
#include <signal.h>

#include <ars/utils.h>

namespace ars {

    //-----------------------------------------------------
    // GaussianMixtureEstimator
    //-----------------------------------------------------

    GaussianMixtureEstimator::GaussianMixtureEstimator() :
    gaussians_() { //: means_(), covars_(), weights_() {
    }

    GaussianMixtureEstimator::~GaussianMixtureEstimator() {
    }

    void GaussianMixtureEstimator::clearGaussians() {
        gaussians_.clear();
    }

    size_t GaussianMixtureEstimator::size() const {
        //return weights_.size();
        return gaussians_.size();
    }

    const Vec2d& GaussianMixtureEstimator::mean(int i) const {
        //            ARS_ASSERT(0 <= i && i < means_.size());
        //            return means_.at(i);
        ARS_ASSERT(0 <= i && i < gaussians_.size());
        return gaussians_[i].mean;
    }

    const Mat2d& GaussianMixtureEstimator::covariance(int i) const {
        //            ARS_ASSERT(0 <= i && i < covars_.size());
        //            return covars_.at(i);
        ARS_ASSERT(0 <= i && i < gaussians_.size());
        return gaussians_[i].covar;
    }

    double GaussianMixtureEstimator::weight(int i) const {
        //            ARS_ASSERT(0 <= i && i < weights_.size());
        //            return weights_.at(i);
        ARS_ASSERT(0 <= i && i < gaussians_.size());
        return gaussians_[i].weight;
    }

    const GaussianMixtureEstimator::VectorGaussian& GaussianMixtureEstimator::gaussians() const {
        return gaussians_;
    }

    void GaussianMixtureEstimator::exportGaussians(VecVec2d &means,
            VecMat2d &covariances, std::vector<double> &weights) const {
        int n = gaussians_.size();
        means.resize(n);
        covariances.resize(n);
        weights.resize(n);
        for (int i = 0; i < n; ++i) {
            means[i] = gaussians_[i].mean;
            covariances[i] = gaussians_[i].covar;
            weights[i] = gaussians_[i].weight;
        }
    }

    void GaussianMixtureEstimator::executeEM(const VecVec2d &samples,
            int stepNum) {
        Eigen::MatrixXd c(samples.size(), gaussians_.size());
        //Eigen::VectorXd sumC(gaussians_.size());
        double sumC;

        for (int step = 0; step < stepNum; ++step) {

            // Expectation
            for (int i = 0; i < samples.size(); ++i) {
                sumC = 0.0;
                for (int j = 0; j < gaussians_.size(); ++j) {
                    c(i, j) = gaussians_[j].weight * gaussians_[j].eval(samples[i]);
                    sumC += c(i, j);
                }
                for (int j = 0; j < gaussians_.size(); ++j) {
                    c(i, j) = c(i, j) / sumC;
                }
            }

            // Maximization
            for (int j = 0; j < gaussians_.size(); ++j) {
                gaussians_[j].mean = Vec2d::Zero();
                gaussians_[j].covar = Mat2d::Zero();
                gaussians_[j].weight = 0.0;
                sumC = 0.0;

                // Computes the mean value
                for (int i = 0; i < samples.size(); ++i) {
                    gaussians_[j].mean += samples[i] * c(i, j);
                    sumC += c(i, j);
                }
                gaussians_[j].mean = gaussians_[j].mean / sumC;
                gaussians_[j].weight = sumC / samples.size();

                // Computes the covariance
                for (int i = 0; i < samples.size(); ++i) {
                    gaussians_[j].covar += (samples[i] - gaussians_[j].mean)
                            * (samples[i] - gaussians_[j].mean).transpose()
                            * c(i, j);
                }
                gaussians_[j].covar = gaussians_[j].covar / sumC;
            }
        }
    }

    //-----------------------------------------------------
    // GaussianMixtureEstimatorScan
    //-----------------------------------------------------

    GaussianMixtureEstimatorScan::GaussianMixtureEstimatorScan() :
    GaussianMixtureEstimator(), intervals_(), distanceGap_(0.8), distanceSplit_(
    0.5), sigmaMin_(0.1) {
    }

    GaussianMixtureEstimatorScan::~GaussianMixtureEstimatorScan() {
    }

    void GaussianMixtureEstimatorScan::compute(const VecVec2d &samples) {
        Vec2d mean;
        Mat2d covar;
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
                //ARS_PRINT("interv [" << interv.first << ", " << interv.last << "] num " << interv.num);
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
                //ARS_PRINT("interv [" << interv.first << ", " << interv.last << "] num " << interv.num);
                Gaussian g;
                //estimateGaussianFromPoints(samples, interv.first, interv.last, g.mean, g.covar);
                estimateGaussianFromSegment(samples, interv.first, interv.last,
                        g.mean, g.covar);
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

    void GaussianMixtureEstimatorScan::findFarthest(const VecVec2d &points,
            int first, int last, int &farthest, double &distMax) const {
        Vec2d dp;
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

    void GaussianMixtureEstimatorScan::estimateGaussianFromPoints(
            const VecVec2d &points, int first, int last, Vec2d &mean,
            Mat2d &covar) const {
        Mat2d l, v;
        Vec2d tmp;
        double sigmaMinSquare = sigmaMin_ * sigmaMin_;

        ARS_ASSERT(first >= 0 && last < points.size());

        // Computes the mean value vector
        mean = Vec2d::Zero();
        for (int i = first; i <= last; ++i) {
            mean += points[i];
        }
        mean = mean / (last - first + 1);

        // Computes the covariance
        covar = Mat2d::Zero();
        for (int i = first; i <= last; ++i) {
            tmp = (points[i] - mean);
            covar = tmp * tmp.transpose();
        }

        if (first == last) {
            // Only one point: use the point uncertainty
            covar << sigmaMinSquare, 0.0, 0.0, sigmaMinSquare;
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

    void GaussianMixtureEstimatorScan::estimateGaussianFromSegment(
            const VecVec2d &points, int first, int last, Vec2d &mean,
            Mat2d &covar) const {
        Mat2d v;
        Vec2d tmp;
        double lmin, lmax, theta, dfirst, dlast;
        double sigmaMinSquare = sigmaMin_ * sigmaMin_;

        ARS_ASSERT(first >= 0 && last < points.size());

        // Computes the mean value vector
        mean = Vec2d::Zero();
        for (int i = first; i <= last; ++i) {
            mean += points[i];
        }
        mean = mean / (last - first + 1);

        // Computes the covariance
        covar = Mat2d::Zero();
        for (int i = first; i <= last; ++i) {
            tmp = (points[i] - mean);
            covar = tmp * tmp.transpose();
        }

        if (first == last) {
            // Only one point: use the point uncertainty
            covar << sigmaMinSquare, 0.0, 0.0, sigmaMinSquare;
        } else {
            covar = covar / (last - first);
            diagonalize(covar, lmin, lmax, theta);
            dfirst = cos(theta) * points[first](0) + sin(theta) * points[first](1);
            dlast = cos(theta) * points[last](0) + sin(theta) * points[last](1);
            lmax = 0.2 * (dlast - dfirst) * (dlast - dfirst);
            if (lmin < sigmaMinSquare) {
                lmin = sigmaMinSquare;
            }
            if (lmax < sigmaMinSquare) {
                lmax = sigmaMinSquare;
            }
            covar << lmax, 0.0, 0.0, lmin;
            v = Eigen::Rotation2Dd(theta);
            covar = v * covar * v.transpose();
            //            ARS_PRINT("[" << first << "," << last << "]: lmin " << lmin << ", lmax " << lmax << ", theta[deg] " << (180.0 / M_PI * theta)
            //                    << "\ncovar\n" << covar);
            //            diagonalize(covar, l, v);
            //            if (l(0, 0) < sigmaMinSquare)
            //                l(0, 0) = sigmaMinSquare;
            //            if (l(1, 1) < sigmaMinSquare)
            //                l(1, 1) = sigmaMinSquare;
            //            covar = v * l * v.transpose();
        }
    }

    //-----------------------------------------------------
    // GaussianMixtureEstimatorMeanShift
    //-----------------------------------------------------

    GaussianMixtureEstimatorMeanShift::GaussianMixtureEstimatorMeanShift() :
    kernelNum_(0), sigmaMin_(0.1), clusterDist_(4.0), meanShiftTol_(2.0), iterationNumMax_(
    30) {
    }

    GaussianMixtureEstimatorMeanShift::~GaussianMixtureEstimatorMeanShift() {
    }

    void GaussianMixtureEstimatorMeanShift::compute(const VecVec2d &samples) {
        VecVec2d meansCurr = samples;
        VecVec2d meansNext(samples.size());
        DisjointSet clusterLabels;
        std::vector<double> clusterIntraDistanceMax;
        std::vector<DisjointSet::id_type> clusterIds;
        Gaussian gaussianKernel;
        double intraDistMaxMax, sigmaMinSquare;
        int iterNum;

        sigmaMinSquare = sigmaMin_ * sigmaMin_;

        ARS_ASSERT(meanShiftTol_ < clusterDist_);

        intraDistMaxMax = std::numeric_limits<double>::max();
        iterNum = 0;
        //	ARS_VARIABLE4(intraDistMaxMax, meanShiftTol_, iterNum, iterationNumMax_);

        while (intraDistMaxMax > meanShiftTol_ && iterNum < iterationNumMax_) {
            updateMeans(meansCurr, meansNext, clusterLabels,
                    clusterIntraDistanceMax);

            clusterIds.clear();
            clusterLabels.parents(std::back_inserter(clusterIds));
            ARS_PRINT(
                    "Iteration " << iterNum << ": " << clusterIds.size() << " clusters: ");
            intraDistMaxMax = 0.0;
            for (auto &cid : clusterIds) {
                //                std::cout << "  cluster " << cid << " samples " << clusterLabels.size(cid) << ": "
                //                        << "max intra-distance " << clusterIntraDistanceMax[cid] << "\n";
                if (clusterIntraDistanceMax[cid] > intraDistMaxMax) {
                    intraDistMaxMax = clusterIntraDistanceMax[cid];
                }
            }
            std::swap(meansCurr, meansNext);
            iterNum++;
        }
        ARS_PRINT("Found " << clusterLabels.setNum() << " clusters");

        // Computes the parameters of the Gaussian kernels associated to each cluster.
        // To find the index of the Gaussian, unfortunately there is a complex index mapping.
        //    samples[i] -> clusterId = clusterLabels.find(i) -> position of clusterId in clusterIds[]
        gaussians_.resize(clusterIds.size());
        for (auto &g : gaussians_) {
            g.mean = Vec2d::Zero();
            g.covar = Mat2d::Zero();
            g.weight = 0.0;
        }
        for (int i = 0; i < samples.size(); ++i) {
            int sampleId = clusterLabels.find(i);
            auto it = std::lower_bound(std::begin(clusterIds), std::end(clusterIds),
                    sampleId);
            size_t gaussianIndex = std::distance(std::begin(clusterIds), it);
            if (it != std::end(clusterIds) && 0 <= gaussianIndex
                    && gaussianIndex < gaussians_.size()) {
                gaussians_[gaussianIndex].mean += samples[i];
                gaussians_[gaussianIndex].covar += samples[i]
                        * samples[i].transpose();
                gaussians_[gaussianIndex].weight += 1.0;
            }
        }
        for (auto &g : gaussians_) {
            if (g.weight > 1.0) {
                g.mean = g.mean / g.weight;
                g.covar = (g.covar - g.weight * g.mean * g.mean.transpose())
                        / (g.weight - 1.0);
                g.weight = g.weight / samples.size();
                saturateEigenvalues(g.covar, sigmaMinSquare);
            } else {
                g.mean = g.mean / g.weight;
                g.covar << sigmaMinSquare, 0.0, 0.0, sigmaMinSquare;
                g.weight = g.weight / samples.size();
            }
        }

        //executeEM(samples, iterationNumMax_);
    }

    void GaussianMixtureEstimatorMeanShift::updateMeans(
            const VecVec2d &meansCurr, VecVec2d &meansNext,
            DisjointSet &clusterLabels,
            std::vector<double> &clusterIntraDistMax) const {
        int n = meansCurr.size();
        std::vector<double> weights(n, 1.0); // evalKernel(0.0) = 1.0
        double d, w;

        clusterLabels.initSets(n);
        clusterIntraDistMax.resize(n);
        std::fill(std::begin(clusterIntraDistMax), std::end(clusterIntraDistMax),
                0.0);

        // Update rule of means:
        //
        //   meansNext[i] = (\sum_{j} meansCurr[j] * w(i,j)) / (\sum_{j} w(i,j))
        //
        // where w(i,j) = K(meansCurr[j] - meansCurr[i]) and w(i,j) = w(j,i).
        // To avoid to recompute the kernel of pairwise distance twice,
        // when it evaluates pair (i,j) both the i-th and the j-th terms are updated.
        //    weights[i] += w(i,j);
        for (int i = 0; i < n; ++i) {
            meansNext[i] = meansCurr[i];
            for (int j = i + 1; j < n; ++j) {
                d = (meansCurr[j] - meansCurr[i]).norm() / (2.0 * sigmaMin_);
                w = exp(-d * d);
                meansNext[i] += meansCurr[j] * w;
                weights[i] += w;
                meansNext[j] += meansCurr[i] * w;
                weights[j] += w;

                // UPDATE of clusters and intra-cluster max distance.
                // a) Case items i and j belong to the same cluster.
                //    It updates the maximum intra-cluster distance.
                // b) Case join clusters of items i and j if their distance is
                //    less than threshold.
                int clusterI = clusterLabels.find(i);
                int clusterJ = clusterLabels.find(j);
                if (clusterI == clusterJ) {
                    if (d > clusterIntraDistMax[clusterI]) {
                        clusterIntraDistMax[clusterI] = d;
                    }
                } else if (d < clusterDist_) {
                    double distMaxI = clusterIntraDistMax[clusterI];
                    double distMaxJ = clusterIntraDistMax[clusterJ];
                    int clusterJoined = clusterLabels.join(i, j);
                    if (d > distMaxI && d > distMaxJ) {
                        clusterIntraDistMax[clusterJoined] = d;
                    } else if (distMaxI > distMaxJ) {
                        clusterIntraDistMax[clusterJoined] = distMaxI;
                    } else {
                        clusterIntraDistMax[clusterJoined] = distMaxJ;
                    }
                }
            }
            meansNext[i] = meansNext[i] / weights[i];
        }
    }

    //-----------------------------------------------------
    // GaussianMixtureEstimatorHierarchical
    //-----------------------------------------------------

    GaussianMixtureEstimatorHierarchical::GaussianMixtureEstimatorHierarchical() :
    data_(), sigmaMin_(0.1), covarWidth_(0.2), chi2Thres_(5.99146), iseThres_(0.3), inlierPerc_(0.60), levelMax_(32) {
    }

    GaussianMixtureEstimatorHierarchical::~GaussianMixtureEstimatorHierarchical() {
    }

    double GaussianMixtureEstimatorHierarchical::getSigmaMin() const {
        return sigmaMin_;
    }

    double GaussianMixtureEstimatorHierarchical::getIseThreshold() const {
        return iseThres_;
    }

    void GaussianMixtureEstimatorHierarchical::setSigmaMin(double sm) {
        sigmaMin_ = sm;
        data_.setRes(sigmaMin_);
    }

    void GaussianMixtureEstimatorHierarchical::setCovarWidth(double cw) {
        covarWidth_ = cw;
    }

    void GaussianMixtureEstimatorHierarchical::setInlierPerc(double ip) {
        inlierPerc_ = ip;
        if (inlierPerc_ < 0.0)
            inlierPerc_ = 0.0;
        else if (inlierPerc_ >= 1.0)
            inlierPerc_ = 1.0;
    }

    void GaussianMixtureEstimatorHierarchical::setChiConfidence(double conf) {
        static const int CHI2_DOF = 2;
        ARS_ASSERT(0.0 <= conf && conf <= 1.0);
        boost::math::chi_squared mydist(CHI2_DOF);
        chi2Thres_ = boost::math::quantile(mydist, conf);
    }

    void GaussianMixtureEstimatorHierarchical::setIseThreshold(double iseTh) {
        iseThres_ = iseTh;
    }

    void GaussianMixtureEstimatorHierarchical::setCellSizeMax(double s) {
        levelMax_ = log2Mod((int) ceil(s / sigmaMin_));
    }

    void GaussianMixtureEstimatorHierarchical::compute(
            const VecVec2d &samples) {
        std::deque<ConstInterval> intervals;
        ConstInterval intervalCur;
        ConstIterator mid;
        Gaussian g;
        double num;
        int level;

        num = (double) samples.size();
        if (num == 0) {
            return;
        }

        gaussians_.clear();
        data_.insert(samples);
        intervals.push_back(std::make_pair(std::begin(data_), std::end(data_)));
        while (!intervals.empty()) {
            intervalCur = intervals.front();
            intervals.pop_front();
            level = data_.computeLevel(intervalCur.first, intervalCur.second);
            if (intervalCur.second != intervalCur.first) {
                //			ARS_VARIABLE4(intervalCur.first->index.transpose(),
                //					(intervalCur.second - 1)->index.transpose(), level,
                //					levelMax_);
            }
            if (level <= levelMax_ && estimateGaussianISE(intervalCur.first, intervalCur.second, g.mean, g.covar, g.weight)) {
                gaussians_.push_front(g);
            } else {
                mid = data_.findSplit(intervalCur.first, intervalCur.second);
                //			ARS_PRINT(
                //					"splitting into\n" << "  ([" << intervalCur.first->index.transpose() << "], [" << mid->index.transpose() << "])\n" << "  ([" << mid->index.transpose() << "], [" << (intervalCur.second-1)->index.transpose() << "])");
                if (mid != intervalCur.first && mid != intervalCur.second) {
                    intervals.push_back(std::make_pair(intervalCur.first, mid));
                    intervals.push_back(std::make_pair(mid, intervalCur.second));
                }
            }
        }

    }

    bool GaussianMixtureEstimatorHierarchical::estimateGaussianFromPoints(
            const ConstIterator &beg, const ConstIterator &end, Vec2d &mean,
            Mat2d &covar, double &w) const {
        Mat2d l, v, infoMat;
        Vec2d tmp;
        double sigmaMinSquare = sigmaMin_ * sigmaMin_;
        double distSqr;
        int num, inlier;

        // Computes the mean value vector
        mean = Vec2d::Zero();
        num = 0;
        for (auto it = beg; it != end; ++it) {
            mean += it->value;
            num++;
        }
        mean = mean / num;

        //	ARS_VARIABLE2(num, mean.transpose());

        // Computes the covariance
        covar = Mat2d::Zero();
        for (auto it = beg; it != end; ++it) {
            tmp = (it->value - mean);
            covar += tmp * tmp.transpose();
        }

        if (num <= 1) {
            // Only one point: use the point uncertainty
            covar << sigmaMinSquare, 0.0, 0.0, sigmaMinSquare;
        } else {
            covar = covar / (num - 1);
            diagonalize(covar, l, v);
            if (l(0, 0) < sigmaMinSquare)
                l(0, 0) = sigmaMinSquare;
            if (l(1, 1) < sigmaMinSquare)
                l(1, 1) = sigmaMinSquare;
            covar = v * l * v.transpose();
        }

        inlier = 0;
        infoMat = covar.inverse();
        //ARS_PRINT("covar\n" << covar << "\neigenvalues:\n" << l.transpose() << "\ninfoMat\n" << infoMat);
        for (auto it = beg; it != end; ++it) {
            tmp = (it->value - mean);
            distSqr = tmp.transpose() * infoMat * tmp;
            if (distSqr < chi2Thres_) {
                inlier++;
            }
        }
        //	ARS_VARIABLE2(inlier, chi2Thres_);
        w = num / data_.size();
        return (inlier >= inlierPerc_ * num);
    }

    bool GaussianMixtureEstimatorHierarchical::estimateGaussianFromSegment(
            const ConstIterator &beg, const ConstIterator &end, Vec2d &mean,
            Mat2d &covar, double &w) const {
        Mat2d l, v, infoMat;
        Vec2d tmp;
        double sigmaMinSquare = sigmaMin_ * sigmaMin_;
        int num, inlier;
        double lmin, lmax, theta, ct, st, distSqr, d, dfirst, dlast;

        // Computes the mean value vector
        mean = Vec2d::Zero();
        num = 0;
        for (auto it = beg; it != end; ++it) {
            mean += it->value;
            num++;
        }
        mean = mean / num;

        // Computes the covariance
        covar = Mat2d::Zero();
        for (auto it = beg; it != end; ++it) {
            tmp = (it->value - mean);
            covar += tmp * tmp.transpose();
        }

        if (num <= 1) {
            // Only one point: use the point uncertainty
            covar << sigmaMinSquare, 0.0, 0.0, sigmaMinSquare;
        } else {
            covar = covar / (num - 1);
            diagonalize(covar, lmin, lmax, theta);
            ct = cos(theta);
            st = sin(theta);
            dfirst = 1e+6;
            dlast = -1e+6;
            for (auto it = beg; it != end; ++it) {
                d = ct * it->value(0) + st * it->value(1);
                if (d < dfirst)
                    dfirst = d;
                if (d > dlast)
                    dlast = d;
            }
            lmax = covarWidth_ * (dlast - dfirst) * (dlast - dfirst);
            if (lmin < sigmaMinSquare) {
                lmin = sigmaMinSquare;
            }
            if (lmax < sigmaMinSquare) {
                lmax = sigmaMinSquare;
            }
            covar << lmax, 0.0, 0.0, lmin;
            v = Eigen::Rotation2Dd(theta);
            covar = v * covar * v.transpose();
        }

        inlier = 0;
        infoMat = covar.inverse();
        //ARS_PRINT("covar\n" << covar << "\neigenvalues:\n" << l.transpose() << "\ninfoMat\n" << infoMat);
        for (auto it = beg; it != end; ++it) {
            tmp = (it->value - mean);
            distSqr = tmp.transpose() * infoMat * tmp;
            if (distSqr < chi2Thres_) {
                inlier++;
            }
        }
        //	ARS_VARIABLE2(inlier, chi2Thres_);
        w = num / data_.size();

        return (inlier >= inlierPerc_ * num);
    }

    bool GaussianMixtureEstimatorHierarchical::estimateGaussianISE(const ConstIterator &beg, const ConstIterator &end, Vec2d &mean,
            Mat2d &covar, double &wMerged) const {

        ARS_ASSERT(beg != end);


        double sigmaMinSqrd = sigmaMin_ * sigmaMin_;
        double wOrig = 1.0 / data_.size(); //w_i


        int num = 0;
        mean = Vec2d::Zero();
        for (auto it = beg; it != end; ++it) {
            mean += it->value;
            num++;
        }
        mean /= num;


        //std::cout << "num " << num << std::endl;
        //std::cout << "mean " << mean << std::endl;



        wMerged = wOrig*num; //w_L

        covar = Mat2d::Zero();
        for (auto it = beg; it != end; ++it) {
            Vec2d tmp = it->value - mean;
            covar += tmp * tmp.transpose();
        }
        covar = covar / num;
        covar(0, 0) += sigmaMinSqrd;
        covar(1, 1) += sigmaMinSqrd;

        //std::cout << "covar\n" << covar << std::endl;


        double jNN = 0.0, jLL = 0.0, jNL = 0.0;

        jLL = wMerged * wMerged / (2 * M_PI * sqrt((2 * covar).determinant()));

        double normOrig = 1.0 / (4 * M_PI * sigmaMinSqrd); //normNN
        double gaussConstOrig = -0.25 / sigmaMinSqrd; //const2

        Mat2d covarNL = covar + sigmaMinSqrd * Mat2d::Identity();
        Mat2d infoNL = covarNL.inverse();
        double normNL = 1.0 / (2 * M_PI * sqrt(covarNL.determinant()));

        for (auto it = beg; it != end; ++it) {
            jNN += 1.0;
            for (auto it2 = it + 1; it2 != end; ++it2) {
                jNN += 2 * exp(gaussConstOrig * (it->value - it2->value).dot(it->value - it2->value));
            }

            jNL += exp(-0.5 * (it->value - mean).transpose() * infoNL * (it->value - mean));
        }
        jNN *= wOrig * wOrig * normOrig;

        jNL *= wOrig * wMerged * normNL;


        double nise = 1 - 2 * jNL / (jNN + jLL);

        if (num == 1) {
            //std::cout << "nise " << nise << std::endl;
            //            std::cout << "wMerged " << wMerged << std::endl;
            //                        std::cout << "wOrig " << wOrig << std::endl;
            //std::cout << "mean\n" << mean.transpose() << std::endl;
            //std::cout << "beg " << (beg->value).transpose() << std::endl;
            //std::cout << "normOrig " << normOrig << std::endl;
            //std::cout << "normNL " << normNL << std::endl;
            //std::cout << "jNN " << jNN << std::endl;
            //std::cout << "jNL " << jNL << std::endl;
            //std::cout << "jLL " << jLL << std::endl;
            return true;
            //            std::cin.get();
        }

        //        std::cout << "--------------\n\n\n" << std::endl;

        return nise < iseThres_;
    }



}// end of namespace

