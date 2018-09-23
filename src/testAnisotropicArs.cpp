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
#include <iostream>
#include <Eigen/Dense>

//#include <ars/ArsKernel.h>
#include <ars/definitions.h>

void getMatrix(double lmin, double lmax, double ang, ars::Matrix2& mat);

int main(int argc, char** argv) {
    ars::Vector2 mean1, mean2;
    ars::Matrix2 covar1, covar2;

    mean1 << 1.0, 3.0;
    mean2 << 4.1, 1.5;
    getMatrix(1.0, 2.0, M_PI / 180.0 * 30.0, covar1);
    getMatrix(1.0, 2.0, M_PI / 180.0 * (-40.0), covar2);

    std::cout << "mean1: " << mean1.transpose() << "\ncovar1\n" << covar1 << std::endl;
    std::cout << "mean2: " << mean2.transpose() << "\ncovar2\n" << covar2 << std::endl;

    return 0;
}

void getMatrix(double lmin, double lmax, double ang, ars::Matrix2& mat) {
    Eigen::Rotation2Dd rot(ang);
    mat << lmin, 0.0,
            0.0, lmax;
    mat = rot * mat * rot.inverse();
    std::cout << "creating matrix with eigenvalues: lmin " << lmin << ", lmax " << lmax << "\n"
            << "  angle " << ang << " rad -> " << (180.0 / M_PI * ang) << std::endl;

    Eigen::SelfAdjointEigenSolver<ars::Matrix2> eigensolver(mat);
    double angle1 = fmod(atan2(eigensolver.eigenvectors().col(0)(1), eigensolver.eigenvectors().col(0)(0)) + 2.0 * M_PI, M_PI);
    double angle2 = fmod(atan2(eigensolver.eigenvectors().col(1)(1), eigensolver.eigenvectors().col(1)(0)) + 2.0 * M_PI, M_PI);
    std::cout << "matrix\n" << mat << std::endl
            << "eigenvalues: " << eigensolver.eigenvalues().transpose()
            << "\neigenvectors\n" << eigensolver.eigenvectors() << "\n"
            << "eigenvector 1 angle " << angle1 << " rad -> " << (180.0 / M_PI * angle1) << " deg\n"
            << "eigenvector 2 angle " << angle2 << " rad -> " << (180.0 / M_PI * angle2) << " deg\n"
            << std::endl;
}

