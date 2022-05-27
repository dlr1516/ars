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

    //    TranslationRefiner::TranslationRefiner() : pointsSrc_(std::vector<Vector2>), pointsDst_(std::vector<Vector2>) {
    //    }

    TranslationRefiner::TranslationRefiner(VectorVector2& ptsSrc, VectorVector2& ptsDst, const Eigen::Translation2d& transl) :
    pointsSrc_(ptsSrc), pointsDst_(ptsDst) {

        Eigen::Affine2d transf;
        transf.setIdentity();
        transf.translation() = transl.translation();

        lastTransf_ = transf;
        for (Vector2& pt : pointsSrc_)
            pt += transl.translation();

    }

    TranslationRefiner::TranslationRefiner(VectorVector2& ptsSrc, VectorVector2& ptsDst, const Eigen::Affine2d& transf) :
    pointsSrc_(ptsSrc), pointsDst_(ptsDst) {

        lastTransf_ = transf;
        for (Vector2& pt : pointsSrc_)
            pt = transf * pt;
    }

    TranslationRefiner::~TranslationRefiner() {
    }

    int TranslationRefiner::computeAssociations() {
        if (associations_.size() > 0) //TODO: maybe this if() can be moved elsewhere
            associations_.clear();

        const int ptsSrcSz = pointsSrc_.size();
        const int ptsDstSz = pointsDst_.size();
        int newAssociations = 0;
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
            if (correspIdx != -1) {
                newAssociations++;
                IndicesPair p(i, correspIdx);
                associations_.push_back(p); //TODO: pushback is good just for first time computing associations!! Need to implement the else{} of associate()
            }
        }
        return newAssociations;
    }

    bool TranslationRefiner::associate(int iteration) {
        std::cout << "Translation Refiner: associate" << std::endl;

        if (iteration == 0 || associations_.empty()) {
            numNewAssocLast_ = computeAssociations();
            return true;
        } else {
            numNewAssocLast_ = computeAssociations();
            if (numNewAssocLast_ / associations_.size() < minNewAssocPerc_) //non-varying associations stopping condition
                return false;
            return true;
        }
    }

    void TranslationRefiner::computeProcrustes(Eigen::Affine2d &transf) {
        //TODO: revise these (specifically associated points handling) steps

        // Procustes
        ars::Vector2 meanSrc, meanDst;
        Eigen::MatrixXd S = Eigen::MatrixXd::Zero(2, 2);
        int n = std::min(pointsSrc_.size(), pointsDst_.size()); // works for now, counting on dummy examples not to have fake correspondences
        int numAssociations = 0;
        int sizeDst = pointsDst_.size();
        meanSrc = ars::Vector2::Zero(2, 1);
        meanDst = ars::Vector2::Zero(2, 1);
        for (int i = 0; i < n; ++i) {
            // if (pointsDst_[i].associated) {
            meanSrc += pointsSrc_[i];
            meanDst += pointsDst_[i];
            numAssociations++;
            // }
        }
        //        std::cout << "numAssociations " << numAssociations << std::endl;
        if (numAssociations != 0) {
            meanSrc = meanSrc * 1.0f / numAssociations;
            meanDst = meanDst * 1.0f / numAssociations;
        } else {
            return;
        }

        for (int i = 0; i < n && i < sizeDst; ++i) {
            // if (!pointsDst_[i].associated)
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

        ars::Vector2 transl = ars::Vector2::Zero(2, 1);
        transl += meanDst - R * meanSrc;

        transf.linear() = R;
        transf.translation() = (meanDst - R * meanSrc);
        transf.makeAffine();
    }

    void TranslationRefiner::icp(Eigen::Affine2d& transfOut) {
        std::cout << "Translation Refiner: ICP" << std::endl;
        double dist = std::numeric_limits<double>::max();
        //TODO: transf pointsSrc points (appropriately) before each iteration of ICP
        for (int i = 0; (i < maxIterations_) && (dist < stopThresh_); ++i) {
            if (associate(i)) {
                computeProcrustes(lastTransf_);
            } else
                break;
        }
        transfOut = lastTransf_;
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