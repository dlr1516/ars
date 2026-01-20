/**
 * ARS - Angular Radon Spectrum
 * Copyright (C) 2025 Dario Lodi Rizzini - Ernesto Fontana.
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
 * */
#include <ars/HoughTranslationEstimator.h>

namespace ars{

HoughTranslationEstimator::HoughTranslationEstimator(double rhoStep, double rhoMax, double sigma) 
    : rhoStep_(rhoStep), rhoMax_(rhoMax), sigma_(sigma) {
        
        invRhoStep_ = 1/rhoStep_;
        rhoNum_ = 2*(rhoMax_*invRhoStep_);
        rhoThresh_ = 3*sigma_;
        initLut();
}

void HoughTranslationEstimator::init(double rhoStep, double rhoMax, double sigma){
    rhoStep_ = rhoStep;
    invRhoStep_ = 1/rhoStep_;
    rhoMax_ = rhoMax;
    sigma_ = sigma;
    rhoNum_ = 2*(rhoMax_*invRhoStep_);
    rhoThresh_ = 3*sigma_;
    initLut();
}

void HoughTranslationEstimator::compute(const VectorVector2 &source, const VectorVector2 &target, Translation& xTrans, Translation& yTrans){
    computeX(source, target, xTrans);
    computeY(source, target, yTrans);
}

void HoughTranslationEstimator::initLut()
{
    double sigma2 = sigma_*sigma_;
    double aDen = 1/(std::sqrt(2*M_PI*sigma2));
    double eDen = 1/(2*sigma2);

    int steps = rhoThresh_*invRhoStep_;
    expLut_.resize(steps);
    //since the gaussian is symmetric, only the right side is calculated
    for(int i = 0; i < steps; i++){
        double rho2 = std::pow(i*rhoStep_, 2);
        double e = aDen*std::exp(-(rho2*eDen));
        expLut_[i] = e;
    }
}

void HoughTranslationEstimator::computeX(const VectorVector2 &source, const VectorVector2 &target, Translation& trans){
    houghSrcX_.assign(rhoNum_, .0);
    houghDstX_.assign(rhoNum_, .0);
    for(auto& point: source){
        double x = point[0];
        updateHough(x, houghSrcX_);
    }

    for(auto& point: target){
        double x = point[0];
        updateHough(x, houghDstX_);
    }
    findCorrelationMax(houghSrcX_, houghDstX_, trans);
}

void HoughTranslationEstimator::computeY(const VectorVector2 &source, const VectorVector2 &target, Translation& trans){
    houghSrcY_.assign(rhoNum_, .0);
    houghDstY_.assign(rhoNum_, .0);
    for(auto& point: source){
        double y = point[1];
        updateHough(y, houghSrcY_);
    }

    for(auto& point: target){
        double y = point[1];
        updateHough(y, houghDstY_);
    }
    findCorrelationMax(houghSrcY_, houghDstY_, trans);
}

void HoughTranslationEstimator::updateHough(double rho, std::vector<double>& hough){
    int firstIdx = std::floor((rho+rhoMax_)*invRhoStep_);
    int lutIdx = 0;
    for(int i = 0; i < expLut_.size(); i++){
        int idx = i+firstIdx;
        int idxSym = firstIdx-lutIdx;
        if(idx >= 0){
            if(idx < hough.size())          hough[idx] += expLut_[lutIdx];
            if(lutIdx != 0 && idxSym > 0 && 
                idxSym < hough.size())      hough[idxSym] += expLut_[lutIdx];
        }
        lutIdx++;
    }
}

void HoughTranslationEstimator::findCorrelationMax(std::vector<double> &src, std::vector<double> &dst, Translation& trans){
    double t = .0;
    double maxCorr = -1.0;
    for(int i = 0; i < rhoNum_; i++){
        double corr = .0;
        int startIdx = i-(rhoNum_*.5);
        for(int j = 0; j < dst.size(); j++){  
            int idx = startIdx+j;
            if(idx >= 0 && idx < src.size())    corr += std::sqrt(src[idx]*dst[j]);
        }
        if(corr > maxCorr){
            maxCorr = corr;
            t = (i * rhoStep_) - rhoMax_;
        }
    }

    trans.translation = t;
    trans.correlation = maxCorr;
}

void HoughTranslationEstimator::plotXHistograms(std::ostream& out){
    for (int i = 0; i < houghSrcX_.size(); i++) {
        auto val = houghSrcX_[i];
        double x = (i * rhoStep_) - rhoMax_;
        out << x << " " << val << "\n";
    }
    out << "e" << std::endl;
    for (int i = 0; i < houghDstX_.size(); i++) {
        auto val = houghDstX_[i];
        double x = (i * rhoStep_) - rhoMax_;
        out << x << " " << val << "\n";
    }
    out << "e" << std::endl;
}

void HoughTranslationEstimator::plotYHistograms(std::ostream& out){
        for (int i = 0; i < houghSrcY_.size(); i++) {
            auto val = houghSrcY_[i];
            double x = (i * rhoStep_) - rhoMax_;
            out << x << " " << val << "\n";
        }
        out << "e" << std::endl;
        for (int i = 0; i < houghDstY_.size(); i++) {
            auto val = houghDstY_[i];
            double x = (i * rhoStep_) - rhoMax_;
            out << x << " " << val << "\n";
        }
        out << "e" << std::endl;
}

}//end of namespace
