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
        a = 0.5 * (m.z - m.w);
        b = 0.5 * (m.x + m.y);
        //ARS_VARIABLE2(a, b);

        theta = 0.5 * atan2(-b, a);

        c = cos(theta);
        s = sin(theta);
        lmax = m.w * c * c + m.z * s * s + (m.x + m.y) * c * s;
        lmin = m.w * s * s + m.z * c * c - (m.x + m.y) * c * s;
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
        cuars::resetToZero(l);
        cuars::setDiagonal(l, lmax, lmin);
        cuars::make2dRotMat(v, theta);
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
        fillRowMajor(covar, lmax, 0.0, 0.0, lmin);
        make2dRotMat(v, theta);
        //                covar = v * covar * v.transpose();
        Mat2d tmpProdResult;
        mat2dProd(tmpProdResult, v, covar);
        transpose(v);
        mat2dProd(covar, tmpProdResult, v);
        transpose(v); //transpose back after using it for the product
    }

    int ceilPow2(int n) {
        ARS_ASSERT(n > 0);

        int exponent = ceil(log2(n));

        int nPadded = std::pow<int>(2, exponent);
        std::cout << "ceilPow2(" << n << ") = " << nPadded << std::endl;

        return nPadded;
    }

    int sumNaturalsUpToN(int n) {
        ARS_ASSERT(n > 0);

        int result = 0.5 * n * (n + 1);

        std::cout << "sumNaturals(" << n << ") = " << result << std::endl;

        return result;
    }

    // --------------------------------------------------------
    // Below: Vec2d and Mat2d util functions (simpler reimplementation of basic Eigen functions)
    // --------------------------------------------------------

    void resetToZero(Vec2d& vec) {
        vec.x = 0.0;
        vec.y = 0.0;
    }

    void resetToZero(Mat2d& mtx) {
        mtx.w = 0.0;
        mtx.x = 0.0;
        mtx.y = 0.0;
        mtx.z = 0.0;
    }

    void setToIdentity(Mat2d& mtx) {
        //        data_[0 * Two + 0] = 1.0; // = data[0]
        //        data_[0 * Two + 1] = 0.0; // = data[1]
        //        data_[1 * Two + 0] = 0.0; // = data[2]
        //        data_[1 * Two + 1] = 1.0; // = data[3]
        mtx.w = 1.0;
        mtx.x = 0.0;
        mtx.y = 0.0;
        mtx.z = 1.0;
    }

    void setDiagonal(Mat2d& mtx, double a11, double a22) {
        mtx.w = a11;
        mtx.x = 0.0;
        mtx.y = 0.0;
        mtx.z = a22;
    }

    void make2dRotMat(Mat2d& mtx, double theta) {
        double cth = cos(theta); //avoiding useless function calling
        double sth = sin(theta);
        mtx.w = cth;
        mtx.x = -sth;
        mtx.y = sth;
        mtx.z = cth;
    }

    void fillRowMajor(Mat2d& mtx, double a, double b, double c, double d) {
        mtx.w = a;
        mtx.x = b;
        mtx.y = c;
        mtx.z = d;
    }

    void scalarMul(Vec2d& vec, double d) {
        vec.x *= d;
        vec.y *= d;
    }

    Vec2d scalarMulWRV(const Vec2d& vec, double d) {
        Vec2d res;

        res.x = vec.x * d;
        res.y = vec.y * d;

        return res;
    }

    void scalarMul(Mat2d& mtx, double d) {
        mtx.w *= d;
        mtx.x *= d;
        mtx.y *= d;
        mtx.z *= d;
    }

    Mat2d scalarMulWRV(const Mat2d& mtx, double d) {
        Mat2d res;

        res.w = mtx.w * d;
        res.x = mtx.x * d;
        res.y = mtx.y * d;
        res.z = mtx.z * d;

        return res;
    }

    void scalarDiv(Vec2d& vec, double d) {
        if (d == 0)
            ARS_ASSERT(false);

        vec.x /= d;
        vec.y /= d;
    }

    Vec2d scalarDivWRV(const Vec2d& vec, double d) {
        if (d == 0)
            ARS_ASSERT(false);

        Vec2d res;

        res.x = vec.x / d;
        res.y = vec.y / d;

        return res;
    }

    void scalarDiv(Mat2d& mtx, double d) {
        if (d == 0)
            ARS_ASSERT(false);

        mtx.w /= d;
        mtx.x /= d;
        mtx.y /= d;
        mtx.z /= d;
    }

    Mat2d scalarDivWRV(const Mat2d& mtx, double d) {
        if (d == 0)
            ARS_ASSERT(false);

        Mat2d res;

        res.w = mtx.w / d;
        res.x = mtx.x / d;
        res.y = mtx.y / d;
        res.z = mtx.z / d;

        return res;
    }

    void transpose(Mat2d& mtx) {
        double tmp = mtx.x;
        mtx.x = mtx.y;
        mtx.y = tmp;
    }

    Mat2d transposeWRV(const Mat2d& mtx) {
        Mat2d res;
        res.x = mtx.y;
        res.y = mtx.x;
        return res;
    }

    double mat2dDeterminant(const Mat2d& mtx) {
        //        return data_[0 * Two + 0] * data_[1 * Two + 1] - data_[0 * Two + 1] * data_[1 * Two + 0];
        return mtx.w * mtx.z - mtx.x * mtx.y;
    }

    double mat2dTrace(const Mat2d& mtx) {
        //        return data_[0 * Two + 0] + data_[1 * Two + 1];
        return mtx.w + mtx.z;
    }

    void mat2dInvert(Mat2d& mtx) {
        double detInv = 1.0 / mat2dDeterminant(mtx);

        double aOrig = mtx.w;
        double bOrig = mtx.x;
        double cOrig = mtx.y;
        double dOrig = mtx.z;

        mtx.w = dOrig * detInv;
        mtx.x = -bOrig * detInv;
        mtx.y = -cOrig * detInv;
        mtx.z = aOrig * detInv;
    }

    Mat2d mat2dInverse(const Mat2d& mtx) {
        double detInv = 1.0 / mat2dDeterminant(mtx);

        double aOrig = mtx.w;
        double bOrig = mtx.x;
        double cOrig = mtx.y;
        double dOrig = mtx.z;

        Mat2d r;

        r.w = dOrig * detInv;
        r.x = -bOrig * detInv;
        r.y = -cOrig * detInv;
        r.z = aOrig * detInv;

        return r;
    }

    void mat2dSum(Mat2d& resultMtx, const Mat2d& aMtx, const Mat2d& bMtx) {
        resultMtx.w = aMtx.w + bMtx.w;
        resultMtx.x = aMtx.x + bMtx.x;
        resultMtx.y = aMtx.y + bMtx.y;
        resultMtx.z = aMtx.z * bMtx.z;
    }

    Mat2d mat2dSumWRV(const Mat2d& aMtx, const Mat2d& bMtx) {
        Mat2d resultMtx;
        resultMtx.w = aMtx.w + bMtx.w;
        resultMtx.x = aMtx.x + bMtx.x;
        resultMtx.y = aMtx.y + bMtx.y;
        resultMtx.z = aMtx.z + bMtx.z;
        return resultMtx;
    }

    void mat2dDiff(Mat2d& resultMtx, const Mat2d& aMtx, const Mat2d& bMtx) {
        resultMtx.w = aMtx.w - bMtx.w;
        resultMtx.x = aMtx.x - bMtx.x;
        resultMtx.y = aMtx.y - bMtx.y;
        resultMtx.z = aMtx.z - bMtx.z;
    }

    Mat2d mat2dDiffWRV(const Mat2d& aMtx, const Mat2d& bMtx) {
        Mat2d resultMtx;
        resultMtx.w = aMtx.w - bMtx.w;
        resultMtx.x = aMtx.x - bMtx.x;
        resultMtx.y = aMtx.y - bMtx.y;
        resultMtx.z = aMtx.z - bMtx.z;
        return resultMtx;
    }

    void mat2dPlusEq(Mat2d& resultMtx, const Mat2d& aMtx) {
        resultMtx.w += aMtx.w;
        resultMtx.x += aMtx.x;
        resultMtx.y += aMtx.y;
        resultMtx.z += aMtx.z;
    }

    void mat2dProd(Mat2d& resultMtx, const Mat2d& aMtx, const Mat2d& bMtx) {
        resultMtx.w = aMtx.w * bMtx.w + aMtx.x * bMtx.y;
        resultMtx.x = aMtx.w * bMtx.x + aMtx.x * bMtx.z;
        resultMtx.y = aMtx.y * bMtx.w + aMtx.z * bMtx.y;
        resultMtx.z = aMtx.x * bMtx.y + aMtx.z * bMtx.z;
    }

    Mat2d mat2dProdWRV(const Mat2d& aMtx, const Mat2d& bMtx) {
        Mat2d resultMtx;
        resultMtx.w = aMtx.w * bMtx.w + aMtx.x * bMtx.y;
        resultMtx.x = aMtx.w * bMtx.x + aMtx.x * bMtx.z;
        resultMtx.y = aMtx.y * bMtx.w + aMtx.z * bMtx.y;
        resultMtx.z = aMtx.x * bMtx.y + aMtx.z * bMtx.z;
        return resultMtx;
    }

    void threeMats2dProd(Mat2d& resultMtx, const Mat2d& aMtx, const Mat2d& bMtx, const Mat2d& cMtx) {
        Mat2d tmp;

        //        mat2dProd(tmp, aMtx, bMtx);
        //        mat2dProd(resultMtx, tmp, cMtx);

        Mat2d aCopy = aMtx;
        Mat2d bCopy = bMtx;
        Mat2d cCopy = cMtx;

        mat2dProd(tmp, aCopy, bCopy);
        mat2dProd(resultMtx, tmp, cCopy);
    }

    double vec2norm(const Vec2d& v) {
        return sqrt(v.x * v.x + v.y * v.y);
    }

    void vec2sum(Vec2d& result, const Vec2d& a, const Vec2d& b) {
        result.x = a.x + b.x;
        result.y = a.y + b.y;
    }

    Vec2d vec2sumWRV(const Vec2d& a, const Vec2d& b) {
        Vec2d result;
        result.x = a.x + b.x;
        result.y = a.y + b.y;
        return result;
    }

    void vec2dPlusEq(Vec2d& result, const Vec2d& v) {
        result.x += v.x;
        result.y += v.y;
    }

    void vec2diff(Vec2d& result, const Vec2d& a, const Vec2d& b) {
        result.x = a.x - b.x;
        result.y = a.y - b.y;
    }

    Vec2d vec2diffWRV(const Vec2d& a, const Vec2d& b) {
        Vec2d result;
        result.x = a.x - b.x;
        result.y = a.y - b.y;
        return result;
    }

    double vec2dotProduct(const Vec2d& a, const Vec2d& b) {
        return a.x * b.x + a.y * b.y;
    }

    void vec2outerProduct(Mat2d& result, const Vec2d& a, const Vec2d& b) {
        result.w = a.x * b.x;
        result.x = a.y * b.x;
        result.y = a.x * b.y;
        result.z = a.y * b.y;
    }

    Mat2d vec2outerProductWRV(const Vec2d& a, const Vec2d& b) {
        Mat2d result;

        result.w = a.x * b.x;
        result.x = a.y * b.x;
        result.y = a.x * b.y;
        result.z = a.y * b.y;

        return result;
    }

    Vec2d row2VecTimesMat2WRV(const Vec2d& v, const Mat2d& m) {
        Vec2d result;
        result.x = (v.x * m.w) + (v.y * m.y);
        result.y = (v.x * m.x) + (v.y * m.z);

        //        result.isCol_ = false;

        return result;
    }

    //affine matrices related

    void preTransfVec2(Vec2d& p, const Affine2d& t) {
        if (t.isLastRowOK()) {
            p.x = (p.x * t.data_[0 * Three + 0]) + (p.y * t.data_[0 * Three + 1]) + (t.data_[0 * Three + 2]);
            p.y = (p.x * t.data_[1 * Three + 0]) + (p.y * t.data_[1 * Three + 1]) + (t.data_[1 * Three + 2]);
            //p.z = 1.0;
        } else {
            printf("ERROR: Transf Matrix affine scale != 1\n");
        }
    }

    void preRotateAff2(Affine2d& t, double angle) {
        if (t.isLastRowOK()) {
            double cth = cos(angle);
            double sth = sin(angle);
            Affine2d tTmpCopy = t;
            //first row
            t.data_[0 * cuars::Three + 0] = (tTmpCopy.data_[0 * cuars::Three + 0] * cth) - (tTmpCopy.data_[1 * cuars::Three + 0] * sth);
            t.data_[0 * cuars::Three + 1] = (tTmpCopy.data_[0 * cuars::Three + 1] * cth) - (tTmpCopy.data_[1 * cuars::Three + 1] * sth);
            t.data_[0 * cuars::Three + 2] = (tTmpCopy.data_[0 * cuars::Three + 2] * cth) - (tTmpCopy.data_[1 * cuars::Three + 2] * sth);
            //second row
            t.data_[1 * cuars::Three + 0] = (tTmpCopy.data_[0 * cuars::Three + 0] * sth) + (tTmpCopy.data_[1 * cuars::Three + 0] * cth);
            t.data_[1 * cuars::Three + 1] = (tTmpCopy.data_[0 * cuars::Three + 1] * sth) + (tTmpCopy.data_[1 * cuars::Three + 1] * cth);
            t.data_[1 * cuars::Three + 2] = (tTmpCopy.data_[0 * cuars::Three + 2] * sth) + (tTmpCopy.data_[1 * cuars::Three + 2] * cth);

            // third (last) row should already be ok
            //            t.data_[2 * cuars::Three + 0] = 0.0;
            //            t.data_[2 * cuars::Three + 1] = 0.0;
            //            t.data_[2 * cuars::Three + 2] = 1.0;

            //            std::cout << "t after prerotation" << std::endl;
            //            std::cout << t;
        } else {
            printf("ERROR: Transf Matrix affine scale != 1\n");
        }
    }

    void preTranslateAff2(Affine2d& t, double x, double y) {
        if (t.isLastRowOK()) {
            //just last column: the other two remain untouched
            t.data_[0 * cuars::Three + 2] = t.data_[0 * cuars::Three + 2] + x; // += x
            t.data_[1 * cuars::Three + 2] = t.data_[1 * cuars::Three + 2] + y; // += y
            //            t.data_[2 * cuars::Three + 2] = 1.0;
        } else {
            printf("ERROR: Transf Matrix affine scale != 1\n");
        }
    }

    void preTranslateAff2(Affine2d& t, const Vec2d& p) {
        if (t.isLastRowOK()) {
            //just last column: the other two remain untouched
            t.data_[0 * cuars::Three + 2] = t.data_[0 * cuars::Three + 2] + p.x; // += x
            t.data_[1 * cuars::Three + 2] = t.data_[1 * cuars::Three + 2] + p.y; // += y
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
            out.data_[0 * cuars::Three + 2] = (a.at(0, 0) * b.at(0, 2)) + (a.at(0, 1) * b.at(1, 2)) + a.at(0, 2); // + a.at(0,2) * b.at(2,2) = a.at(0,2) because b.at(2,2) = 1.0
            out.data_[1 * cuars::Three + 2] = (a.at(1, 0) * b.at(0, 2)) + (a.at(1, 1) * b.at(1, 2)) + a.at(1, 2); // + a.at(1,2) * b.at(2,2) = a.at(1,2) because b.at(2,2) = 1.0
            //            out.data_[2 * cuars::Three + 2] = (a.at(2, 0) * b.at(0, 2)) + (a.at(2, 1) * b.at(1, 2)) + (a.at(2, 2) + b.at(2, 2));
            out.data_[2 * cuars::Three + 2] = 1.0;

        } else {
            printf("ERROR: Transf Matrix last row != 0  0  1\n");
        }
        return out;
    }

    Vec2d aff2TimesVec2WRV(const Affine2d& mAff, const Vec2d& p) {
        Vec2d result;
        if (mAff.isLastRowOK()) {
            result.x = (mAff.at(0, 0) * p.x) + (mAff.at(0, 1) * p.y) + (mAff.at(0, 2));
            result.y = (mAff.at(1, 0) * p.x) + (mAff.at(1, 1) * p.y) + (mAff.at(1, 2));
            //result scale factor = 1
        } else {
            printf("ERROR: Transf Matrix last row != 0  0  1\n");
        }
        return result;
    }

} // end of namespace

