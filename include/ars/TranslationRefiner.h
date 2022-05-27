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

namespace ars {

    //    template <size_t Dim, typename Scalar = double>

    class TranslationRefiner {
    public:
        ars::Vector2 transl_;
        //ars::VectorVector2 pointsSrcTransl_;
        //ars::VectorVector2 pointsDst_;


        /**
         * Default constructor
         */
        //        TranslationRefiner();

        /**
         * Constructor thought for use with ARS rot + transl estimation: translates of @param transl each point in pointsSrc point set
         * that supposedly has been already rotated before using ConsensusTranslationEstimator
         */
        TranslationRefiner(VectorVector2& ptsSrc, VectorVector2& ptsDst, const Eigen::Translation2d& transl);

        /**
         * More generic constructor, that applies a full transformation according to @param transf
         * to ptsSrc
         */
        TranslationRefiner(VectorVector2& ptsSrc, VectorVector2& ptsDst, const Eigen::Affine2d& transf);

        /**
         * Default constructor
         */
        virtual ~TranslationRefiner();

        /**
         * Iterates association and Procrustes problem solving until stopping condition is reached
         * 2 stopping conditions are implemented: 
         *  - max number of iterations
         *  - computed transformation changes below a certain threshold from one iteration to the next
         */
        void icp(Eigen::Affine2d& transfOut);

        /**
         * Setter for maxIterations_ private member
         * Needed for stopping condition
         */
        void setMaxIterations(int numMaxIter);

        /**
         * Setter for stopThresh_ private member
         * Needed for stopping condition
         */
        void setStopThresh(double th);

        /**
         * Setter for association distance threshold: associations are considered only when distance between points
         * is less than @param aDistTh
         */
        void setAssocDistTh(double aDistTh = 0.1);

        /**
         * Setter for the minimum percentage of associations that has to be changed from previous associating iteration
         * is less than @param aDistTh
         */
        void setMinNewAssocPerc(double perc = 0.1);


    private:
        /**
         * Method used inside associate() to compute/re-compute the associations
         * Returns the number of new associations
         */
        int computeAssociations();


        /**
         * Associate each translated point, according to initial/current guess, to closest point in dst
         * @return True when proceeding with Procrusted (associations have varied substantially enough); False otherwise
         */
        bool associate(int iteration);

        /**
         * Solve Procrustes problem finding the most appropriate matrix that links @param pointsSrc and @param pointsDst
         * point sets, based on associations computed at previous step.
         * Saves the result affine matrix in reference @param transf
         */
        void computeProcrustes(Eigen::Affine2d& transf);

        /*
         * Private Members
         */
        Eigen::Affine2d lastTransf_;

        int maxIterations_;
        double stopThresh_;

        double assocDistTh_;

        VectorVector2 &pointsSrc_;
        VectorVector2 &pointsDst_;

        int numNewAssocLast_;
        using IndicesPair = std::pair<int, int>;
        using IndicesPairVec = std::vector<IndicesPair>;
        IndicesPairVec associations_;
        double minNewAssocPerc_;
    };
}

#endif /* ARS_TRANSLATION_REFINER_H */

