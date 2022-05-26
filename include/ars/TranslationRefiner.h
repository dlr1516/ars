/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.h to edit this template
 */

/* 
 * File:   TranslationRefiner.h
 * Author: rimlab
 *
 * Created on May 26, 2022, 5:06 PM
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
        TranslationRefiner();

        /**
         * Default constructor
         */
        virtual ~TranslationRefiner();

        /**
         * Associate each translated point, according to initial/current guess, to closest point in dst
         */
        void associate();

        /**
         * Solve Procrustes problem finding the most appropriate matrix that links @param pointsSrc and @param pointsDst
         * point sets, based on associations computed at previous step.
         * Saves the result affine matrix in reference @param transf
         */
        void computeProcrustes(const ars::VectorVector2 &pointsSrc, const ars::VectorVector2 &pointsDst, Eigen::Affine2d &transf);

        /**
         * Iterates association and Procrustes problem solving until stopping condition is reached
         * 2 stopping conditions are implemented: 
         *  - max number of iterations
         *  - computed transformation changes below a certain threshold from one iteration to the next
         */
        void icp();

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

    private:
        int maxIterations_;
        double stopThresh_;

    };
}

#endif /* ARS_TRANSLATION_REFINER_H */

