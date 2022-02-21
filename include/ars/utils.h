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
#ifndef UTILS_H
#define UTILS_H

#include <ars/definitions.h>

#include <cmath>

namespace cuars {

    /**
     * Computes the diagonalization of the given positive definite matrix m. 
     * 
     *   m = rot(theta) * diag([lmax lmin]) rot(-theta)
     * 
     * @param m the input positive definite matrix
     * @param lmin the minimum eigenvalue
     * @param lmax the maximum eigenvalue
     * @param theta the angle of the eigenvector corresponding to lmax w.r.t. axis x
     */
    void diagonalize(const Mat2d& m, double& lmin, double& lmax, double& theta);

    /**
     * Computes the diagonalization of the given positive definite matrix m. 
     * The relation among the matrices:
     * 
     *   m = v * l * v.transpose()
     * 
     * @param m the input positive definite matrix
     * @param l the matrix of eigenvalues (the maximum eigenvalue first)
     * @param v the matrix with eigenvectors on columns
     */
    void diagonalize(const Mat2d& m, Mat2d& l, Mat2d& v);

    /**
     * Saturates the eigenvalues of the input covariance matrix. 
     * @param covar
     * @param sigmaMinSquare
     */
    void saturateEigenvalues(Mat2d& covar, double sigmaMinSquare);

    // Below: Vec2d and Mat2d util functions (simpler reimplementation of basic Eigen functions)

    void resetToZero(Vec2d& vec);

    void resetToZero(Mat2d& mtx);

    void setDiagonal(Mat2d& mtx, double a11, double a22);

    void make2dRotMat(Mat2d& mtx, double theta);

    void fillRowMajor(Mat2d& mtx, double a, double b, double c, double d);

    void transpose(Mat2d& mtx);

    void mat2dSum(Mat2d& resultMtx, const Mat2d& aMtx, const Mat2d& bMtx);

    void mat2dPlusEq(Mat2d& resultMtx, const Mat2d& aMtx);

    void mat2dProd(Mat2d& resultMtx, const Mat2d& aMtx, const Mat2d & bMtx);

    void threeMats2dProd(Mat2d& resultMtx, const Mat2d& aMtx, const Mat2d& bMtx, const Mat2d & cMtx);

    double vec2norm(const Vec2d & v);

    void vec2sum(Vec2d& result, const Vec2d& a, const Vec2d & b);

    Vec2d vec2sumWRV(const Vec2d& a, const Vec2d & b);

    void vec2diff(Vec2d& result, const Vec2d& a, const Vec2d & b);

    Vec2d vec2diffWRV(const Vec2d& a, const Vec2d & b);

    double vec2dotProduct(Vec2d& a, Vec2d & b);

    void vec2outerProduct(Mat2d& result, const Vec2d& a, const Vec2d & b); //"anti-dot" product: terms are switched

    Mat2d vec2outerProductWRV(Vec2d& a, Vec2d & b); //"anti-dot" product: terms are switched

    Vec2d row2VecTimesMat2WRV(Vec2d& v, Mat2d & m);

    //affine matrices related

    void preTransfVec2(Vec2d& p, const Affine2d & t);

    void preRotateAff2(Affine2d& t, double angle);

    void preTranslateAff2(Affine2d& t, double x, double y);

    void preTranslateAff2(Affine2d& t, const Vec2d& p);

    void aff2Prod(Affine2d& out, const Affine2d& a, const Affine2d & b);

    Affine2d aff2ProdWRV(const Affine2d& a, const Affine2d & b);

    Vec2d aff2TimesVec2WRV(const Affine2d& mAff, const Vec2d& p);




} // end of namespace

#endif /* UTILS_H */

