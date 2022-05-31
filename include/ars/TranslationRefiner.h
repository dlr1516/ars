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

#ifndef ARS_TRANSLATION_REFINER_H
#define ARS_TRANSLATION_REFINER_H

#include <ars/definitions.h>
#include <Eigen/Dense>

namespace ars
{

    //    template <size_t Dim, typename Scalar = double>

    class TranslationRefiner
    {
    public:
        /**
         * Default constructor. Not really useful because of no good way to initialize point vector references...
         */
        TranslationRefiner();

        /**
         * Constructor thought for use with ARS rot + transl estimation: translates of @param transl each point in pointsSrc point set
         * that supposedly has been already rotated before using ConsensusTranslationEstimator
         */
        TranslationRefiner(VectorVector2 &ptsSrc, VectorVector2 &ptsDst, const Eigen::Translation2d &transl);

        /**
         * Analogous of above constructor with transl, just this time with rotation only
         */
        TranslationRefiner(VectorVector2 &ptsSrc, VectorVector2 &ptsDst, const Eigen::Rotation2Dd &rot);

        /**
         * More generic constructor, that applies a full transformation according to @param transf
         * to ptsSrc
         */
        TranslationRefiner(VectorVector2 &ptsSrc, VectorVector2 &ptsDst, const Eigen::Affine2d &transf);

        /**
         * Default constructor
         */
        virtual ~TranslationRefiner();

        /**
         * Iterates association and Procrustes problem solving until stopping condition is reached
         * 3 stopping conditions are implemented:
         *  - max number of iterations
         *  - computed transformation changes below a certain threshold from one iteration to the next
         *  - associations from one iteration to the next don't vary substantially
         */
        void icp(Eigen::Affine2d &transfOut);

        /**
         * Simple solving of Procrustes problem on 2 point sets where pointsSrc_[i] is given as associated
         * to pointsDst_[i], and the size of the two is the same
         */
        void icpNoAssoc(Eigen::Affine2d &transfOut);

        /**
         * Setter for maxIterations_ private member
         * Needed for stopping condition
         */
        void setMaxIterations(int numMaxIter);

        /**
         * Setter for stopMatDistTh_ private member
         * Needed for stopping condition
         */
        void setStopMatDistTh(double th);

        /**
         * Setter for the minimum percentage of associations that has to be changed from previous associating iteration
         * is less than @param aDistTh
         */
        void setMinNewAssocRatio(double ratio);

    private:
        /**
         * Transform each point in pointsSrc set according to transf_ member
         */
        void transformPtsSrc();

        /**
         * Transform each point in pointsSrc set according to @param transf member
         * after setting transf_ member equal to @param transf
         */
        void transformPtsSrcAfterProcrustes(const Eigen::Affine2d &transf);

        /**
         * Method used inside associate() to compute the associations for the first time
         * Returns the number of good associations found
         */
        int computeAssociations();

        /**
         * Method used inside associate() to update the associations in steps of ICP after the first
         * Returns the number of new associations
         */
        int updateAssociations();

        /**
         * Associate each translated point, according to initial/current guess, to closest point in dst
         * @return True when proceeding with Procrusted (associations have varied substantially enough); False otherwise
         */
        bool associate(int iteration);

        /**
         * Solve Procrustes problem finding the most appropriate matrix that links @param pointsSrc and @param pointsDst
         * point sets, based on associations computed at previous step.
         * Saves the result affine matrix in function member transf
         */
        void computeProcrustes();

        /**
         * Setter for association distance threshold: associations are considered only when distance between points
         * is less than @param aDistTh
         */
        void setAssocDistTh(double aDistTh = 10);

        /*
         * Private Members
         */
        Eigen::Affine2d transf_; // serves also as initial guess

        //stopping conditions thresholds
        int maxIterations_;
        double stopMatDistTh_;
        double minNewAssocRatio_;

        //associated pts cannot have distance > than assocDistTh_
        double assocDistTh_;

        VectorVector2 &pointsSrc_;
        VectorVector2 &pointsDst_;


        int numNewAssoc_;
        int numRealAssoc_;
        using IndicesPair = std::pair<int, int>;
        using IndicesPairVec = std::vector<IndicesPair>;
        IndicesPairVec associations_;
    };
}

#endif /* ARS_TRANSLATION_REFINER_H */
