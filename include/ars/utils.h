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

namespace ars {
    
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
    void diagonalize(const Matrix2& m, double& lmin, double& lmax, double& theta);
    
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
    void diagonalize(const Matrix2& m, Matrix2& l, Matrix2& v);
    
    /**
     * Saturates the eigenvalues of the input covariance matrix. 
     * @param covar
     * @param sigmaMinSquare
     */
    void saturateEigenvalues(Matrix2& covar, double sigmaMinSquare);
    
} // end of namespace

#endif /* UTILS_H */

