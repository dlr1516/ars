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

namespace ars {

    void diagonalize(const Matrix2& m, double& lmin, double& lmax, double& theta) {
        double a, b, c, s;

        // Diagonalizes sigma12
        a = 0.5 * (m(1, 1) - m(0, 0));
        b = 0.5 * (m(0, 1) + m(1, 0));
        //ARS_VARIABLE2(a, b);

        theta = 0.5 * atan2(-b, a);

        c = cos(theta);
        s = sin(theta);
        lmax = m(0, 0) * c * c + m(1, 1) * s * s + (m(0, 1) + m(1, 0)) * c * s;
        lmin = m(0, 0) * s * s + m(1, 1) * c * c - (m(0, 1) + m(1, 0)) * c * s;
        //ARS_VARIABLE3(theta, lmax, lmin);

        if (lmax < lmin) {
            theta += 0.5 * M_PI;
            std::swap(lmax, lmin);
            //ARS_PRINT("after swap: lmin " << lmin << " < lmax " << lmax << ", theta + PI/2: " << theta);
        }
    }

    void diagonalize(const Matrix2& m, Matrix2& l, Matrix2& v) {
        double lmin, lmax, theta;

        diagonalize(m, lmin, lmax, theta);
        l = Matrix2::Zero();
        l.diagonal() << lmax, lmin;
        v = Eigen::Rotation2Dd(theta);
    }

    void saturateEigenvalues(Matrix2& covar, double sigmaMinSquare) {
        Matrix2 v;
        double lmin, lmax, theta;

        diagonalize(covar, lmin, lmax, theta);
        if (lmin < sigmaMinSquare) {
            lmin = sigmaMinSquare;
        }
        if (lmax < sigmaMinSquare) {
            lmax = sigmaMinSquare;
        }
        covar << lmax, 0.0,
                0.0, lmin;
        v = Eigen::Rotation2Dd(theta);
        covar = v * covar * v.transpose();
    }

} // end of namespace

