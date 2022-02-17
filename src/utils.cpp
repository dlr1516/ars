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
    
    void mat2dPlusEq(Mat2d& resultMtx, const Mat2d& aMtx) {
        resultMtx.data_[0] += aMtx.data_[0];
        resultMtx.data_[1] += aMtx.data_[1];
        resultMtx.data_[2] += aMtx.data_[2];
        resultMtx.data_[3] += aMtx.data_[3];
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

    double vec2norm(const Vec2d& v) {
        return sqrt(v.x * v.x + v.y * v.y);
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

    void vec2outerProduct(Mat2d& result, const Vec2d& a, const Vec2d& b) {
        result.data_[0 * Two + 0] = a.data_[0] * b.data_[0];
        result.data_[0 * Two + 1] = a.data_[1] * b.data_[0];
        result.data_[1 * Two + 0] = a.data_[0] * b.data_[1];
        result.data_[1 * Two + 1] = a.data_[1] * b.data_[1];
    }

    Mat2d vec2outerProductWRV(const Vec2d& a, const Vec2d& b) {
        Mat2d result;

        result.data_[0 * Two + 0] = a.data_[0] * b.data_[0];
        result.data_[0 * Two + 1] = a.data_[1] * b.data_[0];
        result.data_[1 * Two + 0] = a.data_[0] * b.data_[1];
        result.data_[1 * Two + 1] = a.data_[1] * b.data_[1];

        return result;
    }

    Vec2d row2VecTimesMat2WRV(const Vec2d& v, const Mat2d& m) {
        Vec2d result;
        result.data_[0] = (v.data_[0] * m.data_[0]) + (v.data_[1] * m.data_[2]);
        result.data_[1] = (v.data_[0] * m.data_[1]) + (v.data_[1] * m.data_[3]);

        result.isCol_ = false;

        return result;
    }

    //affine matrices related

    void preTransfVec2(Vec2d& p, const Affine2d& t) {
        if (t.data_[2 * Three + 2] == 1) {
            p.x = (p.x * t.data_[0 * Three + 0]) + (p.y * t.data_[0 * Three + 1]) + (t.data_[0 * Three + 2]);
            p.y = (p.x * t.data_[1 * Three + 0]) + (p.y * t.data_[1 * Three + 1]) + (t.data_[1 * Three + 2]);
            //p.z = 1.0;
        } else {
            printf("ERROR: Transf Matrix affine scale != 1\n");
        }
    }

    void preRotateAff2(Affine2d& t, double angle) {
        if (t.data_[2 * Three + 2] == 1) {
            double cth = cos(angle);
            double sth = sin(theta);
            //first row
            t.data_[0 * cuars::Three + 0] = (t.data_[0 * cuars::Three + 0] * cth) - (t.data_[1 * cuars::Three + 0] * sth);
            t.data_[0 * cuars::Three + 1] = (t.data_[0 * cuars::Three + 1] * cth) - (t.data_[1 * cuars::Three + 1] * sth);
            t.data_[0 * cuars::Three + 2] = (t.data_[0 * cuars::Three + 2] * cth) - (t.data_[1 * cuars::Three + 2] * sth);
            //second row
            t.data_[1 * cuars::Three + 0] = (t.data_[0 * cuars::Three + 0] * sth) + (t.data_[1 * cuars::Three + 0] * cth);
            t.data_[1 * cuars::Three + 1] = (t.data_[0 * cuars::Three + 1] * sth) + (t.data_[1 * cuars::Three + 1] * cth);
            t.data_[1 * cuars::Three + 2] = (t.data_[0 * cuars::Three + 2] * sth) + (t.data_[1 * cuars::Three + 2] * cth);

            // third (last) row should already be ok
            //            t.data_[2 * cuars::Three + 0] = 0.0;
            //            t.data_[2 * cuars::Three + 1] = 0.0;
            //            t.data_[2 * cuars::Three + 2] = 1.0;
        } else {
            printf("ERROR: Transf Matrix affine scale != 1\n");
        }
    }

    void preTranslateAff2(Affine2d& t, double x, double y) {
        if (t.data_[2 * Three + 2] == 1) {
            //just last column: the other two remain untouched
            t.data_[0 * cuars::Three + 2] = t.data_[0 * cuars::Three + 2] + x; // += x
            t.data_[1 * cuars::Three + 2] = t.data_[1 * cuars::Three + 2] + y; // += y
            //            t.data_[2 * cuars::Three + 2] = 1.0;
        } else {
            printf("ERROR: Transf Matrix affine scale != 1\n");
        }
    }

    void aff2Prod(Affine2d& out, const Affine2d& a, const Affine2d& b) {
        if (a.isLastRowOK() && b.isLastRowOK()) {
            //elements not mentioned are implicitly ok (or invariant if they are sum terms, because last rows are [0  0  1])

            //first column
            out.data_[0 * cuars::Three + 0] = (a.at(0, 0) * b.at(0, 0)) + (a.at(0, 1) * b.at(1, 0)); // + a.at(0,2) * b.at(2,0)
            out.data_[1 * cuars::Three + 0] = (a.at(1, 0) * b.at(0, 0)) + (a.at(1, 1) * b.at(1, 0)); // + a.at(1,2) * b.at(2,0)
            //            out.data_[2 * cuars::Three + 0] = (a.at(2, 0) * b.at(0, 0)) + (a.at(2, 1) * b.at(1, 0)) + (a.at(2,2) * b.at(2,0));
            out.data_[2 * cuars::Three + 0] = 0.0;

            //second column
            out.data_[0 * cuars::Three + 1] = (a.at(0, 0) * b.at(0, 1)) + (a.at(0, 1) * b.at(1, 1)); // + a.at(0,2) * b.at(2,1) 
            out.data_[1 * cuars::Three + 1] = (a.at(1, 0) * b.at(0, 1)) + (a.at(1, 1) * b.at(1, 1)); // + a.at(1,2) * b.at(2,1)
            //            out.data_[2 * cuars::Three + 1] = (a.at(2, 0) * b.at(0, 1)) + (a.at(2, 1) * b.at(1, 1)) + (a.at(2,2) * b.at(2,1));
            out.data_[2 * cuars::Three + 1] = 0.0;


            //third column
            out.data_[0 * cuars::Three + 2] = (a.at(0, 0) * b.at(0, 2)) + (a.at(0, 1) * b.at(1, 2)); // + a.at(0,2) * b.at(2,2)
            out.data_[1 * cuars::Three + 2] = (a.at(1, 0) * b.at(0, 2)) + (a.at(1, 1) * b.at(1, 2)); // + a.at(1,2) * b.at(2,2)
            //            out.data_[2 * cuars::Three + 2] = (a.at(2, 0) * b.at(0, 2)) + (a.at(2, 1) * b.at(1, 2)) + (a.at(2, 2) + b.at(2, 2));
            out.data_[2 * cuars::Three + 2] = 1.0;

        } else {
            printf("ERROR: Transf Matrix last row != 0  0  1\n");
        }
    }
    
    void aff2Prod(Affine2d& out, const Affine2d& a, const Affine2d& b) {
        if (a.isLastRowOK() && b.isLastRowOK()) {
            //elements not mentioned are implicitly ok (or invariant if they are sum terms, because last rows are [0  0  1])

            //first column
            out.data_[0 * cuars::Three + 0] = (a.at(0, 0) * b.at(0, 0)) + (a.at(0, 1) * b.at(1, 0)); // + a.at(0,2) * b.at(2,0)
            out.data_[1 * cuars::Three + 0] = (a.at(1, 0) * b.at(0, 0)) + (a.at(1, 1) * b.at(1, 0)); // + a.at(1,2) * b.at(2,0)
            //            out.data_[2 * cuars::Three + 0] = (a.at(2, 0) * b.at(0, 0)) + (a.at(2, 1) * b.at(1, 0)) + (a.at(2,2) * b.at(2,0));
            out.data_[2 * cuars::Three + 0] = 0.0;

            //second column
            out.data_[0 * cuars::Three + 1] = (a.at(0, 0) * b.at(0, 1)) + (a.at(0, 1) * b.at(1, 1)); // + a.at(0,2) * b.at(2,1) 
            out.data_[1 * cuars::Three + 1] = (a.at(1, 0) * b.at(0, 1)) + (a.at(1, 1) * b.at(1, 1)); // + a.at(1,2) * b.at(2,1)
            //            out.data_[2 * cuars::Three + 1] = (a.at(2, 0) * b.at(0, 1)) + (a.at(2, 1) * b.at(1, 1)) + (a.at(2,2) * b.at(2,1));
            out.data_[2 * cuars::Three + 1] = 0.0;


            //third column
            out.data_[0 * cuars::Three + 2] = (a.at(0, 0) * b.at(0, 2)) + (a.at(0, 1) * b.at(1, 2)); // + a.at(0,2) * b.at(2,2)
            out.data_[1 * cuars::Three + 2] = (a.at(1, 0) * b.at(0, 2)) + (a.at(1, 1) * b.at(1, 2)); // + a.at(1,2) * b.at(2,2)
            //            out.data_[2 * cuars::Three + 2] = (a.at(2, 0) * b.at(0, 2)) + (a.at(2, 1) * b.at(1, 2)) + (a.at(2, 2) + b.at(2, 2));
            out.data_[2 * cuars::Three + 2] = 1.0;

        } else {
            printf("ERROR: Transf Matrix last row != 0  0  1\n");
        }
    }
    
    Affine2d aff2ProdWRV(const Affine2d& a, const Affine2d& b) {
        Affine2d out;
        if (a.isLastRowOK() && b.isLastRowOK()) {
            //elements not mentioned are implicitly ok (or invariant if they are sum terms, because last rows are [0  0  1])

            //first column
            out.data_[0 * cuars::Three + 0] = (a.at(0, 0) * b.at(0, 0)) + (a.at(0, 1) * b.at(1, 0)); // + a.at(0,2) * b.at(2,0)
            out.data_[1 * cuars::Three + 0] = (a.at(1, 0) * b.at(0, 0)) + (a.at(1, 1) * b.at(1, 0)); // + a.at(1,2) * b.at(2,0)
            //            out.data_[2 * cuars::Three + 0] = (a.at(2, 0) * b.at(0, 0)) + (a.at(2, 1) * b.at(1, 0)) + (a.at(2,2) * b.at(2,0));
            out.data_[2 * cuars::Three + 0] = 0.0;

            //second column
            out.data_[0 * cuars::Three + 1] = (a.at(0, 0) * b.at(0, 1)) + (a.at(0, 1) * b.at(1, 1)); // + a.at(0,2) * b.at(2,1) 
            out.data_[1 * cuars::Three + 1] = (a.at(1, 0) * b.at(0, 1)) + (a.at(1, 1) * b.at(1, 1)); // + a.at(1,2) * b.at(2,1)
            //            out.data_[2 * cuars::Three + 1] = (a.at(2, 0) * b.at(0, 1)) + (a.at(2, 1) * b.at(1, 1)) + (a.at(2,2) * b.at(2,1));
            out.data_[2 * cuars::Three + 1] = 0.0;


            //third column
            out.data_[0 * cuars::Three + 2] = (a.at(0, 0) * b.at(0, 2)) + (a.at(0, 1) * b.at(1, 2)); // + a.at(0,2) * b.at(2,2)
            out.data_[1 * cuars::Three + 2] = (a.at(1, 0) * b.at(0, 2)) + (a.at(1, 1) * b.at(1, 2)); // + a.at(1,2) * b.at(2,2)
            //            out.data_[2 * cuars::Three + 2] = (a.at(2, 0) * b.at(0, 2)) + (a.at(2, 1) * b.at(1, 2)) + (a.at(2, 2) + b.at(2, 2));
            out.data_[2 * cuars::Three + 2] = 1.0;

        } else {
            printf("ERROR: Transf Matrix last row != 0  0  1\n");
        }
        return out;
    }



} // end of namespace

