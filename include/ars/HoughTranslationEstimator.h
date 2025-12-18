/**
 * ARS - Angular Radon Spectrum 
 * Copyright (C) 2017 Dario Lodi Rizzini.
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
#ifndef HOUGH_TRANSLATION_ESTIMATOR_H
#define HOUGH_TRANSLATION_ESTIMATOR_H

#include <iostream>
#include <Eigen/Dense>
#include <ars/ars2d.h>

namespace ars {

    /** 
     * Class for translation estimator between two scan by using the hough transform
     * at angles 0° and 90°, cosidering they are already aligned.
     */
    class HoughTranslationEstimator {
    public:
        struct Translation{
            double translation;
            double correlation;
        };
        /** Constructor with deafult parameters. 
         */
        HoughTranslationEstimator() : expLut_() {};

        /**
         * Constructor with initialization parameters
         * @param rhoStep the dimension of range bin
         * @param rhoMax the maximum value of polar range of lines
         * @param sigma the standard deviation of the sensor used
         */
        HoughTranslationEstimator(double rhoStep, double rhoMax, double sigma);

        /** Default destructor. 
         */
        virtual ~HoughTranslationEstimator() {};

        /**
         * Inits params of the class. 
         * @param rhoStep the dimension of range bin
         * @param rhoMax the maximum value of polar range of lines
         * @param sigma the standard deviation of the sensor used
         */
        void init(double rhoStep, double rhoMax, double sigma);

        /** Computes the estimated translation. 
         */
        void compute(const VectorVector2& source, const VectorVector2& target, Translation& xTrans, Translation& yTrans);

        void plotXHistograms(std::ostream& out);
        void plotYHistograms(std::ostream& out);

    private:
        int rhoNum_;
        double rhoStep_, invRhoStep_;
        double rhoMax_;
        double sigma_;
        double rhoThresh_;
        std::vector<double> houghSrcX_;
        std::vector<double> houghDstX_;
        std::vector<double> houghSrcY_;
        std::vector<double> houghDstY_;

        //clouds are represented as a mixutre of gaussians with the same variance.
        //A lut is used to get the value to use in the transform.
        std::vector<double> expLut_;

        void initLut();

        void computeX(const VectorVector2& source, const VectorVector2& target, Translation& trans);
        void computeY(const VectorVector2& source, const VectorVector2& target, Translation& trans);

        void updateHough(double rho, std::vector<double>& hough);
        void findCorrelationMax(std::vector<double>& src, std::vector<double>& dst, Translation& trans);
    };//end of class
} //end of namespace
#endif //HOUGH_TRANSLATION_ESTIMATOR_H