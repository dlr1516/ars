/**
 * ARS - Angular Radon Spectrum 
 * Copyright (C) 2017-2020 Dario Lodi Rizzini
 *               2021- Dario Lodi Rizzini, Ernesto Fontana
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

#include <ars/TranslationRefiner.h>

namespace ars {

    //    TranslationRefiner::TranslationRefiner() : pointsSrc_(), pointsDst_(), transl(Eigen::Affine2d::Identity()) {
    //    }

    TranslationRefiner::TranslationRefiner(VectorVector2& ptsSrc, VectorVector2& ptsDst, const Eigen::Translation2d& transl) :
    pointsSrc_(ptsSrc), pointsDst_(ptsDst) {

        Eigen::Affine2d transf(transl);

        transf_ = transf;

        //moved pointsSrc_ update accoring to transf_ initial guess in transformPtsSrc() method
        //        for (Vector2& pt : pointsSrc_)
        //            pt += transl.translation();
    }

    TranslationRefiner::TranslationRefiner(VectorVector2& ptsSrc, VectorVector2& ptsDst, const Eigen::Rotation2Dd& rot) :
    pointsSrc_(ptsSrc), pointsDst_(ptsDst) {

        Eigen::Affine2d transf(rot);

        transf_ = transf;

        //moved pointsSrc_ update accoring to transf_ initial guess in transformPtsSrc() method
        //        for (Vector2& pt : pointsSrc_)
        //            pt = rot * pt;
    }

    TranslationRefiner::TranslationRefiner(VectorVector2& ptsSrc, VectorVector2& ptsDst, const Eigen::Affine2d& transf) :
    pointsSrc_(ptsSrc), pointsDst_(ptsDst) {

        transf_ = transf;

        //moved pointsSrc_ update accoring to transf_ initial guess in transformPtsSrc() method
        //        for (Vector2& pt : pointsSrc_)
        //            pt = transf * pt;
    }

    TranslationRefiner::~TranslationRefiner() {
    }

    void TranslationRefiner::transformPtsSrc() {
        for (Vector2& pt : pointsSrc_)
            pt = transf_ * pt;
    }

    int TranslationRefiner::computeAssociations() {
        if (associations_.size() > 0) //TODO: maybe this if() can be moved elsewhere (or even just removed)
            associations_.clear();

        const int ptsSrcSz = pointsSrc_.size();
        const int ptsDstSz = pointsDst_.size();
        int goodAssociations = 0;
        for (int i = 0; i < ptsSrcSz; ++i) {
            //            bool associated = false;
            double distBest = std::numeric_limits<double>::max();
            int correspIdx = -1;
            for (int j = 0; j < ptsDstSz; ++j) {
                //TODO maybe: compute centroids directly here

                double distancePtToPt = (pointsDst_.at(j) - pointsSrc_.at(i)).squaredNorm();
                if (distancePtToPt >= assocDistTh_)
                    continue;
                if (distancePtToPt < distBest) {
                    //                    associated = true;
                    distBest = distancePtToPt;
                    correspIdx = j;
                }
            }

            IndicesPair p(i, correspIdx);
            associations_.push_back(p);

            if (correspIdx != -1) { //establishing associations vector "real" size (size without counting bad associations)
                goodAssociations++;
            }
        }
        return goodAssociations;
    }

    int TranslationRefiner::updateAssociations() {
        int newAssociations = 0;

        const int currAssocSz = associations_.size();
        for (int i = 0; i < currAssocSz; ++i) {
            const int idxSrc = associations_[i].first;
            const int idxDst = associations_[i].second; //TODO maybe: if idxDst == -1 -> go next straight away

            const int ptsDstSz = pointsDst_.size();
            double distCurr = (pointsDst_.at(idxDst) - pointsSrc_.at(idxSrc)).squaredNorm();
            double distBest = distCurr;

            int newCorrespIdx = -1;

            for (int j = 0; j < ptsDstSz; ++j) {
                //TODO maybe: compute centroids directly here
                double distPtToPt = (pointsDst_.at(j) - pointsSrc_.at(idxSrc)).squaredNorm();

                if (distPtToPt >= assocDistTh_)
                    continue;
                if (distPtToPt < distBest) {
                    //                    associated = true;
                    distBest = distPtToPt;
                    newCorrespIdx = j;
                }
            }

            if (newCorrespIdx != idxDst && newCorrespIdx != -1) {
                associations_[i].second = newCorrespIdx;
                newAssociations++;
            }
        }

        return newAssociations;
    }

    bool TranslationRefiner::associate(int iteration) {
        std::cout << "Translation Refiner: associate" << std::endl;

        if (iteration == 0 || associations_.empty()) {
            numRealAssoc_ = computeAssociations();
            numNewAssocLast_ = numRealAssoc_;
            return true;
        } else {
            numNewAssocLast_ = updateAssociations();
            if (numNewAssocLast_ / numRealAssoc_ < minNewAssocPerc_) //non-varying associations stopping condition
                return false;
            return true;
        }
    }

    void TranslationRefiner::computeProcrustes() {
        Eigen::Affine2d transf;
        //TODO: revise these steps (specifically: associated points handling)

        // Procustes
        Vector2 meanSrc, meanDst;
        Eigen::MatrixXd S = Eigen::MatrixXd::Zero(2, 2);

        const int n = associations_.size();

        meanSrc = Vector2::Zero(2, 1);
        meanDst = Vector2::Zero(2, 1);
        for (int i = 0; i < n; ++i) {
            const int dstIdx = associations_.at(i).second;

            if (dstIdx == -1)
                continue;

            const int srcIdx = associations_.at(i).first;

            const Vector2 ptSrc = pointsSrc_.at(srcIdx);
            const Vector2 ptDst = pointsDst_.at(dstIdx);

            meanSrc += ptSrc;
            meanDst += ptDst;
        }

        if (n != 0) {
            meanSrc = meanSrc * 1.0f / n;
            meanDst = meanDst * 1.0f / n;
        } else {
            return;
        }

        for (int i = 0; i < n; ++i) {
            const int dstIdx = associations_.at(i).second;

            if (dstIdx == -1)
                continue;

            const int srcIdx = associations_.at(i).first;

            const Vector2 ptSrc = pointsSrc_.at(srcIdx);
            const Vector2 ptDst = pointsDst_.at(dstIdx);

            S += (ptSrc - meanSrc) * (ptDst - meanDst).transpose();
        }
        //        std::cout << "S before SVD" << std::endl
        //                << S << std::endl;

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(S, Eigen::ComputeThinU | Eigen::ComputeThinV);
        float d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
        if (d > 0)
            d = 1.0;
        else
            d = -1.0;
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(2, 2);
        I(1, 1) = d;
        Eigen::MatrixXd R = svd.matrixV() * I * svd.matrixU().transpose();

        //        std::cout << "S after SVD" << std::endl
        //                << S << std::endl;
        //        std::cout << "U after SVD" << std::endl
        //                << S << std::endl;
        //        std::cout << "V after SVD" << std::endl
        //                << S << std::endl;

        Vector2 transl = Vector2::Zero(2, 1);
        transl += meanDst - R * meanSrc;

        transf.linear() = R;
        transf.translation() = (meanDst - R * meanSrc);
        transf.makeAffine();

        transf_ = transf;
    }

    void TranslationRefiner::icp(Eigen::Affine2d & transfOut) {
        std::cout << "Translation Refiner: ICP" << std::endl;
        double dist = std::numeric_limits<double>::max();

        for (int i = 0; (i < maxIterations_) && (dist < stopThresh_); ++i) {
            //            if (i != 0) //maybe needed when points passed as input (through constructor) are alerady transformed
            transformPtsSrc();

            if (associate(i)) {
                const Eigen::Affine2d lastTransf = transf_;

                computeProcrustes();

                dist = (lastTransf.matrix() - transf_.matrix()).norm();
                if (dist < stopThresh_)
                    break;
            } else
                break;
        }
        transfOut = transf_;
    }

    void TranslationRefiner::icpNoAssoc(Eigen::Affine2d & transfOut) {
        Eigen::Affine2d transf;

        // Procustes
        Vector2 meanSrc, meanDst;
        Eigen::MatrixXd S = Eigen::MatrixXd::Zero(2, 2);
        int n = std::min(pointsSrc_.size(), pointsDst_.size()); // works for now, counting on dummy examples not to have fake correspondences
        int numAssociations = 0;
        int sizeDst = pointsDst_.size();
        meanSrc = Vector2::Zero(2, 1);
        meanDst = Vector2::Zero(2, 1);
        for (int i = 0; i < n; ++i) {
            // if (pointsDst[i].associated) {
            meanSrc += pointsSrc_[i];
            meanDst += pointsDst_[i];
            numAssociations++;
            // }
        }
        std::cout << "numAssociations " << numAssociations << std::endl;
        if (numAssociations != 0) {
            meanSrc = meanSrc * 1.0f / numAssociations;
            meanDst = meanDst * 1.0f / numAssociations;
        } else {
            return;
        }

        for (int i = 0; i < n && i < sizeDst; ++i) {
            // if (!pointsDst[i].associated)
            // continue;
            S += (pointsSrc_[i] - meanSrc) * (pointsDst_[i] - meanDst).transpose();
        }
        //        std::cout << "S before SVD" << std::endl
        //                << S << std::endl;

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(S, Eigen::ComputeThinU | Eigen::ComputeThinV);
        float d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
        if (d > 0)
            d = 1.0;
        else
            d = -1.0;
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(2, 2);
        I(1, 1) = d;
        Eigen::MatrixXd R = svd.matrixV() * I * svd.matrixU().transpose();

        //        std::cout << "S after SVD" << std::endl
        //                << S << std::endl;
        //        std::cout << "U after SVD" << std::endl
        //                << S << std::endl;
        //        std::cout << "V after SVD" << std::endl
        //                << S << std::endl;

        Vector2 transl = Vector2::Zero(2, 1);
        transl += meanDst - R * meanSrc;

        transf.linear() = R;
        transf.translation() = (meanDst - R * meanSrc);
        transf.makeAffine();

        transf_ = transf;
        transfOut = transf;
    }

    void TranslationRefiner::setMaxIterations(int numMaxIter) {
        maxIterations_ = numMaxIter;
    }

    void TranslationRefiner::setStopThresh(double th) {
        stopThresh_ = th;
    }

    void TranslationRefiner::setAssocDistTh(double aDistTh) {
        assocDistTh_ = aDistTh;
    }

}