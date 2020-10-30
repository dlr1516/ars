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

#include <ars/definitions.h>
#include <ars/ars2d.h>

double acesRanges[] = {50.00, 50.00, 50.00, 5.26, 5.21, 5.06, 5.01, 3.01, 2.94, 2.89, 2.84, 2.74, 2.69, 2.64, 2.59, 2.54, 2.49, 2.49, 2.44, 2.39, 2.34, 2.29, 2.29, 2.29, 2.39, 2.39, 2.49, 2.51, 2.61, 2.66, 2.76, 2.81, 2.96, 3.01, 3.11, 3.26, 3.01, 3.01, 3.01, 3.06, 3.21, 6.86, 6.86, 6.81, 6.76, 6.71, 6.71, 6.66, 6.61, 6.66, 6.56, 6.56, 6.56, 6.46, 6.46, 6.41, 6.46, 6.46, 4.11, 3.96, 3.96, 4.96, 4.86, 5.21, 7.41, 4.61, 5.16, 6.26, 6.26, 6.31, 4.86, 5.01, 5.86, 5.81, 4.21, 4.26, 4.31, 4.41, 4.39, 4.46, 5.31, 5.06, 5.26, 4.96, 6.01, 5.76, 5.61, 5.36, 5.26, 5.01, 4.21, 4.16, 4.01, 3.91, 3.61, 3.21, 3.26, 3.16, 3.06, 3.01, 3.31, 3.21, 3.16, 2.16, 2.19, 2.16, 2.21, 2.11, 2.01, 2.01, 2.06, 2.84, 2.91, 2.91, 3.01, 3.11, 3.21, 3.81, 4.06, 7.11, 7.06, 7.01, 6.96, 6.86, 4.31, 6.76, 6.71, 6.66, 6.61, 5.46, 5.41, 6.46, 6.21, 6.31, 6.51, 7.26, 7.46, 50.00, 2.01, 1.94, 1.94, 1.94, 2.31, 1.86, 1.84, 1.84, 1.81, 1.96, 26.46, 20.76, 2.11, 2.12, 2.17, 2.14, 2.09, 2.09, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.19, 2.19, 2.24, 2.24, 2.24, 2.24, 2.29, 2.29, 2.29, 2.29, 2.29, 2.39, 2.39, 2.39, 2.44};

void rangeToPoint(double* ranges, int num, double angleMin, double angleRes, ars::VectorVector2& points);



//void getMatrix(double lmin, double lmax, double ang, ars::Matrix2& mat);

int main(int argc, char** argv) {
    ars::AngularRadonSpectrum2d ars;
    ars::VectorVector2 acesPoints;

    rangeToPoint(acesRanges, 180, -0.5 * M_PI, M_PI / 180.0 * 1.0, acesPoints);
    acesPoints.push_back(ars::Vector2::Zero());
    std::cout << "Number of input points: " << acesPoints.size() << std::endl;
    //    ars::Vector2 mean1, mean2;
    //    ars::Matrix2 covar1, covar2;
    //
    //    mean1 << 1.0, 3.0;
    //    mean2 << 4.1, 1.5;
    //    getMatrix(1.0, 2.0, M_PI / 180.0 * 30.0, covar1);
    //    getMatrix(1.0, 2.0, M_PI / 180.0 * (-40.0), covar2);
    //
    //    std::cout << "mean1: " << mean1.transpose() << "\ncovar1\n" << covar1 << std::endl;
    //    std::cout << "mean2: " << mean2.transpose() << "\ncovar2\n" << covar2 << std::endl;

    return 0;
}

void rangeToPoint(double* ranges, int num, double angleMin, double angleRes, ars::VectorVector2& points) {
    Eigen::Vector2d p;
    for (int i = 0; i < num; ++i) {
        double a = angleMin + angleRes * i;
        p << ranges[i] * cos(a), ranges[i] * sin(a);
        points.push_back(p);
    }
}

//void getMatrix(double lmin, double lmax, double ang, ars::Matrix2& mat) {
//    Eigen::Rotation2Dd rot(ang);
//    mat << lmin, 0.0,
//            0.0, lmax;
//    mat = rot * mat * rot.inverse();
//    std::cout << "creating matrix with eigenvalues: lmin " << lmin << ", lmax " << lmax << "\n"
//            << "  angle " << ang << " rad -> " << (180.0 / M_PI * ang) << std::endl;
//
//    Eigen::SelfAdjointEigenSolver<ars::Matrix2> eigensolver(mat);
//    double angle1 = fmod(atan2(eigensolver.eigenvectors().col(0)(1), eigensolver.eigenvectors().col(0)(0)) + 2.0 * M_PI, M_PI);
//    double angle2 = fmod(atan2(eigensolver.eigenvectors().col(1)(1), eigensolver.eigenvectors().col(1)(0)) + 2.0 * M_PI, M_PI);
//    std::cout << "matrix\n" << mat << std::endl
//            << "eigenvalues: " << eigensolver.eigenvalues().transpose()
//            << "\neigenvectors\n" << eigensolver.eigenvectors() << "\n"
//            << "eigenvector 1 angle " << angle1 << " rad -> " << (180.0 / M_PI * angle1) << " deg\n"
//            << "eigenvector 2 angle " << angle2 << " rad -> " << (180.0 / M_PI * angle2) << " deg\n"
//            << std::endl;
//}

