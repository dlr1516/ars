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
#include <ars/BBOptimizer1d.h>
#include <ars/HoughSpectrum.h>
#include <ars/HistogramCircularCorrelation.h>


#include <chrono>
#include <ars/thirdparty/gnuplot-iostream.h>

#define PRINT_DIM(X) std::cout << #X << " rows " << X.rows() << " cols " << X.cols() << std::endl;
#define RAD2DEG(X) (180.0/M_PI*(X))

struct BoundInterval {
    double x0;
    double x1;
    double y0;
    double y1;
};

double acesRanges1[] = {50.00, 50.00, 50.00, 5.26, 5.21, 5.06, 5.01, 3.01, 2.94, 2.89, 2.84, 2.74, 2.69, 2.64, 2.59, 2.54, 2.49, 2.49, 2.44, 2.39, 2.34, 2.29, 2.29, 2.29, 2.39, 2.39, 2.49, 2.51, 2.61, 2.66, 2.76, 2.81, 2.96, 3.01, 3.11, 3.26, 3.01, 3.01, 3.01, 3.06, 3.21, 6.86, 6.86, 6.81, 6.76, 6.71, 6.71, 6.66, 6.61, 6.66, 6.56, 6.56, 6.56, 6.46, 6.46, 6.41, 6.46, 6.46, 4.11, 3.96, 3.96, 4.96, 4.86, 5.21, 7.41, 4.61, 5.16, 6.26, 6.26, 6.31, 4.86, 5.01, 5.86, 5.81, 4.21, 4.26, 4.31, 4.41, 4.39, 4.46, 5.31, 5.06, 5.26, 4.96, 6.01, 5.76, 5.61, 5.36, 5.26, 5.01, 4.21, 4.16, 4.01, 3.91, 3.61, 3.21, 3.26, 3.16, 3.06, 3.01, 3.31, 3.21, 3.16, 2.16, 2.19, 2.16, 2.21, 2.11, 2.01, 2.01, 2.06, 2.84, 2.91, 2.91, 3.01, 3.11, 3.21, 3.81, 4.06, 7.11, 7.06, 7.01, 6.96, 6.86, 4.31, 6.76, 6.71, 6.66, 6.61, 5.46, 5.41, 6.46, 6.21, 6.31, 6.51, 7.26, 7.46, 50.00, 2.01, 1.94, 1.94, 1.94, 2.31, 1.86, 1.84, 1.84, 1.81, 1.96, 26.46, 20.76, 2.11, 2.12, 2.17, 2.14, 2.09, 2.09, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.19, 2.19, 2.24, 2.24, 2.24, 2.24, 2.29, 2.29, 2.29, 2.29, 2.29, 2.39, 2.39, 2.39, 2.44};

double acesRanges2[] = {3.04, 2.94, 2.94, 2.84, 2.79, 2.69, 2.64, 2.59, 2.54, 2.49, 2.44, 2.44, 2.39, 2.34, 2.29, 2.24, 2.24, 2.34, 2.34, 2.49, 2.49, 2.56, 2.61, 2.71, 2.81, 2.81, 2.91, 3.06, 3.21, 3.01, 3.01, 3.01, 3.06, 3.16, 3.21, 6.86, 6.81, 6.76, 6.76, 6.71, 6.66, 6.61, 6.61, 6.61, 6.56, 6.56, 6.51, 6.46, 6.41, 4.26, 6.41, 4.16, 4.01, 3.91, 3.81, 4.86, 4.86, 5.11, 4.76, 3.96, 5.46, 6.21, 6.21, 6.26, 4.46, 6.26, 5.76, 5.76, 4.21, 4.26, 4.31, 4.36, 4.39, 4.37, 5.06, 5.06, 6.51, 4.91, 6.01, 5.76, 5.56, 5.36, 5.16, 5.01, 4.21, 4.11, 3.96, 3.86, 3.56, 3.16, 3.21, 3.16, 3.01, 3.01, 3.26, 3.21, 3.11, 2.11, 2.14, 2.11, 2.16, 2.06, 1.96, 1.96, 1.96, 2.71, 2.81, 2.91, 2.96, 3.06, 3.11, 3.76, 3.81, 7.11, 7.06, 7.01, 6.91, 6.81, 6.81, 4.26, 6.71, 6.61, 6.56, 5.81, 5.36, 5.36, 6.16, 6.26, 6.31, 7.31, 7.36, 50.00, 50.00, 1.91, 1.94, 1.89, 1.94, 14.61, 1.81, 1.84, 1.84, 1.79, 1.91, 20.76, 7.11, 2.11, 2.12, 2.12, 2.09, 2.09, 2.09, 2.09, 2.09, 2.09, 2.09, 2.14, 2.09, 2.14, 2.14, 2.19, 2.19, 2.19, 2.19, 2.24, 2.19, 2.24, 2.24, 2.29, 2.24, 2.29, 2.29, 2.34, 2.34, 2.39, 2.44, 2.44, 2.49, 2.49, 2.54, 2.54};


void rangeToPoint(double* ranges, int num, double angleMin, double angleRes, ars::Vector2Vector& points);

void plotBranchBoundBox(std::ostream& out, const std::vector<BoundInterval>& bbbs);

int main() {
    ars::AngularRadonSpectrum2d ars1;
    ars::AngularRadonSpectrum2d ars2;
    ars::HoughSpectrum hs1;
    ars::HoughSpectrum hs2;
    ars::HistogramCircularCorrelation corr;
    ars::Vector2Vector acesPoints1;
    ars::Vector2Vector acesPoints2;
    std::vector<double> correlationFourier;
    std::chrono::system_clock::time_point timeStart, timeStop;
    double sigma = 0.05;
    int fourierOrder = 20;
    double thetaMax1, arsMax1, thetaMax2, arsMax2;
    double thetaTol, fourierTol, thetaMax, corrMax;
    int thnum = 360;
    int hsShiftMax;
    double hsCorrMax;

    ars1.setARSFOrder(fourierOrder);
    ars2.setARSFOrder(fourierOrder);

    rangeToPoint(acesRanges1, 180, -0.5 * M_PI, M_PI / 180.0 * 1.0, acesPoints1);
    acesPoints1.push_back(ars::Vector2::Zero());
    std::cout << "Number of input points in scan 1: " << acesPoints1.size() << std::endl;

    rangeToPoint(acesRanges2, 180, -0.5 * M_PI, M_PI / 180.0 * 1.0, acesPoints2);
    acesPoints2.push_back(ars::Vector2::Zero());
    std::cout << "Number of input points in scan 2: " << acesPoints2.size() << std::endl;


    ars1.initLUT(0.0001);
    ars1.setComputeMode(ars::AngularRadonSpectrum2d::PNEBI_LUT);
    ars1.insertIsotropicGaussians(acesPoints1, sigma);
    arsMax1 = ars1.findMax(thetaMax1);
    std::cout << "ars1: maximum " << arsMax1 << " in " << (180.0 / M_PI * thetaMax1) << " [deg]" << std::endl;

    ars2.initLUT(0.0001);
    ars2.setComputeMode(ars::AngularRadonSpectrum2d::PNEBI_LUT);
    ars2.insertIsotropicGaussians(acesPoints2, sigma);
    arsMax2 = ars2.findMax(thetaMax2);
    std::cout << "ars2: maximum " << arsMax2 << " in " << (180.0 / M_PI * thetaMax2) << " [deg]" << std::endl;

    ars::computeFourierCorr(ars1.coefficients(), ars2.coefficients(), correlationFourier);

    thetaTol = M_PI / 180.0 * 0.2;
    fourierTol = 1.0;
    ars::findGlobalMaxBBFourier(correlationFourier, 0.0, M_PI, thetaTol, fourierTol, thetaMax, corrMax);

    std::cout << "ARS: best correlation for rotation " << (180.0 / M_PI * thetaMax) << " [deg] with max value " << corrMax << std::endl;

    hs1.init(M_PI / thnum, 0.20, 30.0);
    hs2.init(M_PI / thnum, 0.20, 30.0);
    hs1.insertPoint(acesPoints1.begin(), acesPoints1.end());
    hs2.insertPoint(acesPoints2.begin(), acesPoints2.end());
    corr.computeHistSimpleShift(hs1.spectrum(), hs2.spectrum(), thnum, hsShiftMax, hsCorrMax);
    std::cout << "HS: best correlation for rotation " << (180.0 * hsShiftMax / thnum) << " [deg] with max value " << corrMax << std::endl;


    Gnuplot gp("gnuplot -persist");
    //    double vieweps = 5e-3;
    //    //  std::ostream& gp = std::cout;
    gp << "set term wxt 0\n";
    //gp << "set terminal postscript eps enhanced color \"Times-Roman\" 24\n";
    //gp << "set output \"registration_ars.eps\"\n";
    gp << "plot '-' title \"ARS 1\" w l lw 5.0, '-' title \"ARS 2\" w l lw 2.0, '-' title \"ARS corr\" w l lw 5.0\n";
    for (int i = 0; i < thnum; ++i) {
        double a = M_PI / thnum * i;
        gp << (180.0 / thnum * i) << " " << ars1.eval(a) / arsMax1 << "\n";
    }
    gp << "e" << std::endl;
    for (int i = 0; i < thnum; ++i) {
        double a = M_PI / thnum * i;
        gp << (180.0 / thnum * i) << " " << ars2.eval(a) / arsMax2 << "\n";
    }
    gp << "e" << std::endl;
    for (int i = 0; i < thnum; ++i) {
        double a = M_PI / thnum * i;
        gp << (180.0 / thnum * i) << " " << ars::evaluateFourier(correlationFourier, 2.0 * a) / corrMax << "\n";
    }
    gp << "e" << std::endl;

    gp << "set term wxt 1\n";
    //gp << "set output \"registration_ars.eps\"\n";
    gp << "plot '-' title \"HS 1\" w l lw 5.0, '-' title \"HS 2\" w l lw 2.0\n";
    for (int i = 0; i < thnum; ++i) {
        double a = M_PI / thnum * i;
        gp << (180.0 / thnum * i) << " " << hs1.spectrum(a) << "\n";
    }
    gp << "e" << std::endl;
    for (int i = 0; i < thnum; ++i) {
        double a = M_PI / thnum * i;
        gp << (180.0 / thnum * i) << " " << hs2.spectrum(a) << "\n";
    }
    gp << "e" << std::endl;

    gp << "set term wxt 2\n";
    //gp << "set terminal postscript eps enhanced color \"Times-Roman\" 24\n";
    gp << "set output \"registration_input_scans.eps\"\n";
    gp << "set size ratio -1\n";
    gp << "set xrange [0.0:8.0]\n";
    gp << "set yrange [-8.0:8.0]\n";
    gp << "plot '-' title \"scan 1\" w p pt 7 ps 0.5, '-' title \"scan 2\" w p pt 7 ps 0.5\n";
    for (auto& p : acesPoints1) {
        gp << p.x() << " " << p.y() << "\n";
    }
    gp << "e" << std::endl;
    for (auto& p : acesPoints2) {
        gp << p.x() << " " << p.y() << "\n";
    }
    gp << "e" << std::endl;


    return 0;
}

void rangeToPoint(double* ranges, int num, double angleMin, double angleRes, ars::Vector2Vector& points) {
    Eigen::Vector2d p;
    for (int i = 0; i < num; ++i) {
        double a = angleMin + angleRes * i;
        p << ranges[i] * cos(a), ranges[i] * sin(a);
        points.push_back(p);
    }
}

void plotBranchBoundBox(std::ostream& out, const std::vector<BoundInterval>& bbbs) {
    for (auto& bbb : bbbs) {
        out << RAD2DEG(bbb.x0) << " " << bbb.y0 << "\n"
                << RAD2DEG(bbb.x1) << " " << bbb.y0 << "\n"
                << RAD2DEG(bbb.x1) << " " << bbb.y1 << "\n"
                << RAD2DEG(bbb.x0) << " " << bbb.y1 << "\n"
                << RAD2DEG(bbb.x0) << " " << bbb.y0 << "\n\n";
    }
}


