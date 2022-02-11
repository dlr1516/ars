/**
 * ARS - Angular Radon Spectrum 
 * Copyright (C) 2017-2020 Dario Lodi Rizzini.
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
#include <ars/utils.h>
#include <ars/definitions.h>

namespace cuars {

    void diagonalize(const Mat2d& m, double& lmin, double& lmax, double& theta) {
        double a, b, c, s;

        // Diagonalizes sigma12
        a = 0.5 * (m.data_[1 * Two + 1] - m.data_[0 * Two + 0]);
        b = 0.5 * (m.data_[0 * Two + 1] + m.data_[1 * Two + 0]);
        //ARS_VARIABLE2(a, b);

        theta = 0.5 * atan2(-b, a);

        c = cos(theta);
        s = sin(theta);
        lmax = m.data_[0 * Two + 0] * c * c + m.data_[1 * Two + 1] * s * s + (m.data_[0 * Two + 1] + m.data_[1 * Two + 0]) * c * s;
        lmin = m.data_[0 * Two + 0] * s * s + m.data_[1 * Two + 1] * c * c - (m.data_[0 * Two + 1] + m.data_[1 * Two + 0]) * c * s;
        //ARS_VARIABLE3(theta, lmax, lmin);

        if (lmax < lmin) {
            theta += 0.5 * M_PI;
            std::swap(lmax, lmin);
            //ARS_PRINT("after swap: lmin " << lmin << " < lmax " << lmax << ", theta + PI/2: " << theta);
        }
    }

    void diagonalize(const Mat2d& m, Mat2d& l, Mat2d& v) {
        double lmin, lmax, theta;

        diagonalize(m, lmin, lmax, theta);
        l.resetToZero();
        l.setDiagonal(lmax, lmin);
        v.make2dRotMat(theta);
    }

    void saturateEigenvalues(Mat2d& covar, double sigmaMinSquare) {
        Mat2d v;
        double lmin, lmax, theta;

        diagonalize(covar, lmin, lmax, theta);
        if (lmin < sigmaMinSquare) {
            lmin = sigmaMinSquare;
        }
        if (lmax < sigmaMinSquare) {
            lmax = sigmaMinSquare;
        }
        covar.fillRowMajor(lmax, 0.0, 0.0, lmin);
        v.make2dRotMat(theta);
        //                covar = v * covar * v.transpose();
        Mat2d tmpProdResult;
        mat2dProd(tmpProdResult, v, covar);
        v.transpose();
        mat2dProd(covar, tmpProdResult, v);
        v.transpose(); //transpose back after using it for the product
    }

    //

    //    void zeroResetPointerVals(Vec2d& vec) {
    //        vec.data_[0] = 0.0;
    //        vec.data_[1] = 0.0;
    //    }
    //
    //    void zeroResetMatrixVals(Mat2d& mtx) {
    //        mtx.data_[0] = 0.0;
    //        mtx.data_[1] = 0.0;
    //        mtx.data_[2] = 0.0;
    //        mtx.data_[3] = 0.0;
    //    }

    void mat2dSum(Mat2d& resultMtx, const Mat2d& aMtx, const Mat2d& bMtx) {
        resultMtx.data_[0] = aMtx.data_[0] * bMtx.data_[0];
        resultMtx.data_[1] = aMtx.data_[1] * bMtx.data_[1];
        resultMtx.data_[2] = aMtx.data_[2] * bMtx.data_[2];
        resultMtx.data_[3] = aMtx.data_[3] * bMtx.data_[3];
    }

    void mat2dProd(Mat2d& resultMtx, const Mat2d& aMtx, const Mat2d& bMtx) {
        resultMtx.data_[0] = aMtx.data_[0] * bMtx.data_[0] + aMtx.data_[1] * bMtx.data_[2];
        resultMtx.data_[1] = aMtx.data_[0] * bMtx.data_[1] + aMtx.data_[1] * bMtx.data_[3];
        resultMtx.data_[2] = aMtx.data_[2] * bMtx.data_[0] + aMtx.data_[3] * bMtx.data_[2];
        resultMtx.data_[3] = aMtx.data_[1] * bMtx.data_[2] + aMtx.data_[3] * bMtx.data_[3];
    }

    void threeMats2dProd(Mat2d& resultMtx, const Mat2d& aMtx, const Mat2d& bMtx, const Mat2d& cMtx) {
        Mat2d tmp;
        mat2dProd(tmp, aMtx, bMtx);
        mat2dProd(resultMtx, tmp, cMtx);
    }

    void vec2sum(Vec2d& result, const Vec2d& a, const Vec2d& b) {
        result.data_[0] = a.data_[0] + b.data_[0];
        result.data_[1] = a.data_[1] + b.data_[1];
    }

    Vec2d vec2sumWRV(Vec2d& a, Vec2d& b) {
        Vec2d result;
        result.data_[0] = a.data_[0] + b.data_[0];
        result.data_[1] = a.data_[1] + b.data_[1];
        return result;
    }

    void vec2diff(Vec2d& result, const Vec2d& a, const Vec2d& b) {
        result.data_[0] = a.data_[0] - b.data_[0];
        result.data_[1] = a.data_[1] - b.data_[1];
    }

    Vec2d vec2diffWRV(Vec2d& a, Vec2d& b) {
        Vec2d result;
        result.data_[0] = a.data_[0] - b.data_[0];
        result.data_[1] = a.data_[1] - b.data_[1];
        return result;
    }

    double vec2dotProduct(Vec2d& a, Vec2d& b) {
        return a.data_[0] * b.data_[0] + a.data_[1] * b.data_[1];
    }

    void vec2outerProduct(Mat2d& result, Vec2d& a, Vec2d& b) {
        result.data_[0 * Two + 0] = a.data_[0] * b.data_[0];
        result.data_[0 * Two + 1] = a.data_[1] * b.data_[0];
        result.data_[1 * Two + 0] = a.data_[0] * b.data_[1];
        result.data_[1 * Two + 1] = a.data_[1] * b.data_[1];
    }

    Mat2d vec2outerProductWRV(Vec2d& a, Vec2d& b) {
        Mat2d result;

        result.data_[0 * Two + 0] = a.data_[0] * b.data_[0];
        result.data_[0 * Two + 1] = a.data_[1] * b.data_[0];
        result.data_[1 * Two + 0] = a.data_[0] * b.data_[1];
        result.data_[1 * Two + 1] = a.data_[1] * b.data_[1];

        return result;
    }

    Vec2d row2VecTimesMat2WRV(const Vec2d& v, const Mat2d& m) {
        Vec2d result;
        result.data_[0] = v.data_[0] * m.data_[0] + v.data_[1] * m.data_[2];
        result.data_[1] = v.data_[0] * m.data_[1] + v.data_[1] * m.data_[3];

        result.isCol_ = false;

        return result;
    }




} // end of namespace

