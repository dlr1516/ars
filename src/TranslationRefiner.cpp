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

namespace ars
{
    // PUBLIC TranslationRefiner members

    //    TranslationRefiner::TranslationRefiner() : pointsSrc_(), pointsDst_(), transl(Eigen::Affine2d::Identity()) {
    //    }

    TranslationRefiner::TranslationRefiner(VectorVector2 &ptsSrc, VectorVector2 &ptsDst, const Eigen::Translation2d &transl) : pointsSrc_(ptsSrc), pointsDst_(ptsDst)
    {

        Eigen::Affine2d transf(transl);

        transf_ = transf;

        // moved pointsSrc_ update accoring to transf_ initial guess in transformPtsSrc() method
        //         for (Vector2& pt : pointsSrc_)
        //             pt += transl.translation();
    }

    TranslationRefiner::TranslationRefiner(VectorVector2 &ptsSrc, VectorVector2 &ptsDst, const Eigen::Rotation2Dd &rot) : pointsSrc_(ptsSrc), pointsDst_(ptsDst)
    {

        Eigen::Affine2d transf(rot);

        transf_ = transf;

        // moved pointsSrc_ update accoring to transf_ initial guess in transformPtsSrc() method
        //         for (Vector2& pt : pointsSrc_)
        //             pt = rot * pt;
    }

    TranslationRefiner::TranslationRefiner(VectorVector2 &ptsSrc, VectorVector2 &ptsDst, const Eigen::Affine2d &transf) : pointsSrc_(ptsSrc), pointsDst_(ptsDst)
    {

        transf_ = transf;

        // moved pointsSrc_ update accoring to transf_ initial guess in transformPtsSrc() method
        //         for (Vector2& pt : pointsSrc_)
        //             pt = transf * pt;
    }

    TranslationRefiner::~TranslationRefiner()
    {
    }

    void TranslationRefiner::icp(Eigen::Affine2d &transfOut)
    {
        setAssocDistTh();
        std::cout << "Translation Refiner: ICP" << std::endl;
        double dist = std::numeric_limits<double>::max();

        transformPtsSrc(); // transformations are computed incrementally
        for (int i = 0; i < maxIterations_; ++i)
        {
            std::cout << "ICP iteration # " << i << std::endl;
            //            if (i != 0) //maybe needed when points passed as input (through constructor) are alerady transformed

            if (associate(i))
            {
                const Eigen::Affine2d lastTransf = transf_;

                computeProcrustes();

                dist = (lastTransf.matrix() - transf_.matrix()).norm();
                if (dist < stopMatDistTh_)
                {
                    std::cout << "icp: VERY LIMITED CHANGE IN TRANSF -> STOPPING" << std::endl;
                    transfOut = transf_;
                    break;
                }
            }
            else
            {
                std::cout << "icp: ASSOCIATIONS BARELY CHANGE -> STOPPING" << std::endl;
                transfOut = transf_;
                break;
            }

            if (i == maxIterations_ - 1)
            {
                std::cout << "icp: MAX NUM ITERATIONS REACHED -> STOPPING" << std::endl;
                transfOut = transf_;
                break;
            }
        }
    }

    void TranslationRefiner::icpNoAssoc(Eigen::Affine2d &transfOut)
    {
        std::cout << "Translation Refiner: ICP no assoc" << std::endl;

        Eigen::Affine2d transf;

        // Procustes
        Vector2 meanSrc, meanDst;
        Eigen::MatrixXd S = Eigen::MatrixXd::Zero(2, 2);
        int n = std::min(pointsSrc_.size(), pointsDst_.size()); // works for now, counting on dummy examples not to have fake correspondences
        int numAssociations = 0;
        int sizeDst = pointsDst_.size();
        meanSrc = Vector2::Zero(2, 1);
        meanDst = Vector2::Zero(2, 1);
        for (int i = 0; i < n; ++i)
        {
            // if (pointsDst[i].associated) {
            meanSrc += pointsSrc_[i];
            meanDst += pointsDst_[i];
            numAssociations++;
            // }
        }
        std::cout << "numAssociations " << numAssociations << std::endl;
        if (numAssociations != 0)
        {
            meanSrc = meanSrc * 1.0f / numAssociations;
            meanDst = meanDst * 1.0f / numAssociations;
        }
        else
        {
            return;
        }

        for (int i = 0; i < n && i < sizeDst; ++i)
        {
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

    void TranslationRefiner::setMaxIterations(int numMaxIter)
    {
        maxIterations_ = numMaxIter;
    }

    void TranslationRefiner::setStopMatDistTh(double th)
    {
        stopMatDistTh_ = th;
    }

    void TranslationRefiner::setMinNewAssocRatio(double perc)
    {
        minNewAssocRatio_ = perc;
    }

    // PRIVATE TranslationRefiner members

    void TranslationRefiner::transformPtsSrc()
    {
        // std::cout << "applying transf " << std::endl
        //           << transf_.matrix() << std::endl;

        for (Vector2 &pt : pointsSrc_)
            pt = transf_ * pt;
    }

    void TranslationRefiner::transformPtsSrcAfterProcrustes(const Eigen::Affine2d &transf)
    {
        transf_ = transf * transf_;
        transformPtsSrc();
    }

    int TranslationRefiner::computeAssociations()
    {
        if (associations_.size() > 0) // TODO: maybe this if() can be moved elsewhere (or even just removed)
            associations_.clear();

        const int ptsSrcSz = pointsSrc_.size();
        const int ptsDstSz = pointsDst_.size();
        int goodAssociations = 0;
        for (int i = 0; i < ptsSrcSz; ++i)
        {
            //            bool associated = false;
            double distBest = std::numeric_limits<double>::max();
            int correspIdx = -1;
            for (int j = 0; j < ptsDstSz; ++j)
            {
                // TODO maybe: compute centroids directly here

                double distPtToPt = (pointsDst_.at(j) - pointsSrc_.at(i)).squaredNorm();
                if (distPtToPt >= assocDistTh_)
                    continue;
                const double samePtDistTh = 0.001;
                if (distPtToPt < distBest && distPtToPt > samePtDistTh)
                {
                    //                    associated = true;
                    distBest = distPtToPt;
                    correspIdx = j;
                }
            }

            IndicesPair p(i, correspIdx);
            associations_.push_back(p);
            // std::cout << "just computed association " << i << " from point " << pointsSrc_.at(i).transpose() << " to point " << pointsDst_.at(correspIdx).transpose()
            //           << std::endl
            //           << p.first << " " << p.second << std::endl;

            if (correspIdx != -1)
            { // establishing associations vector "real" size (size without counting bad associations)
                goodAssociations++;
            }
        }
        return goodAssociations;
    }

    int TranslationRefiner::updateAssociations()
    {
        int newAssociations = 0;

        const int numAssoc = associations_.size();
        for (int i = 0; i < numAssoc; ++i)
        {
            const int idxSrc = associations_[i].first;
            const int idxDst = associations_[i].second; // TODO maybe: if idxDst == -1 -> go next straight away

            const int ptsDstSz = pointsDst_.size();
            double distCurr = (pointsDst_.at(idxDst) - pointsSrc_.at(idxSrc)).squaredNorm();
            double distBest = distCurr;

            int newCorrespIdx = -1;

            for (int j = 0; j < ptsDstSz; ++j)
            {
                // TODO maybe: compute centroids directly here
                double distPtToPt = (pointsDst_.at(j) - pointsSrc_.at(idxSrc)).squaredNorm();

                if (distPtToPt >= assocDistTh_)
                    continue;
                const double samePtDistTh = 0.001;
                if (distPtToPt < distBest && distPtToPt > samePtDistTh)
                {
                    //                    associated = true;
                    distBest = distPtToPt;
                    newCorrespIdx = j;
                }
            }

            if (newCorrespIdx != idxDst && newCorrespIdx != -1)
            {
                associations_[i].second = newCorrespIdx;
                newAssociations++;
            }
        }

        return newAssociations;
    }

    bool TranslationRefiner::associate(int iteration)
    {
        std::cout << "Translation Refiner: associate" << std::endl;

        if (iteration == 0 || associations_.empty())
        {
            numRealAssoc_ = computeAssociations(); // numRealAssoc_ won't change anymore
            numNewAssoc_ = numRealAssoc_;
            return true;
        }
        else
        {
            numNewAssoc_ = updateAssociations();
            std::cout << "new associations: " << numNewAssoc_ << " out of " << numRealAssoc_ << " total associations" << std::endl;
            const double newAssocRatio = numNewAssoc_ / numRealAssoc_;
            std::cout << "ratio: " << newAssocRatio << std::endl;
            std::cout << "min ratio: " << minNewAssocRatio_ << std::endl;

            if (newAssocRatio < minNewAssocRatio_) // non-varying associations stopping condition
                return false;
            return true;
        }
    }

    void TranslationRefiner::computeProcrustes()
    {
        std::cout << "TranslationRefiner::computeProcrustes()" << std::endl;
        Eigen::Affine2d transf;
        // TODO: revise these steps (specifically: associated points handling)

        // Procustes
        Vector2 meanSrc, meanDst;
        Eigen::MatrixXd S = Eigen::MatrixXd::Zero(2, 2);

        const int n = associations_.size();

        if (n == 0)
            return;

        meanSrc = Vector2::Zero(2, 1);
        meanDst = Vector2::Zero(2, 1);
        for (int i = 0; i < n; ++i)
        {
            const int dstIdx = associations_.at(i).second;

            if (dstIdx == -1)
                continue;

            const int srcIdx = associations_.at(i).first;

            const Vector2 ptSrc = pointsSrc_.at(srcIdx);
            const Vector2 ptDst = pointsDst_.at(dstIdx);

            meanSrc += ptSrc * 1.0f / n;
            meanDst += ptDst * 1.0f / n;
        }

        for (int i = 0; i < n; ++i)
        {
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

        transformPtsSrcAfterProcrustes(transf);

        transf_ = transf * transf_;

        // std::cout << "computed new transf (procrustes output): " << std::endl
        //           << transf_.matrix() << std::endl;
    }

    void TranslationRefiner::setAssocDistTh(double aDistTh)
    {
        assocDistTh_ = aDistTh;
    }

}