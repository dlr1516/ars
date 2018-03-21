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


void rangeToPoint(double* ranges, int num, double angleMin, double angleRes, ars::Point2Vector& points);

void plotBranchBoundBox(std::ostream& out, const std::vector<BoundInterval>& bbbs);

int main() {
    ars::AngularRadonSpectrum2d ars1;
    ars::AngularRadonSpectrum2d ars2;
    ars::Point2Vector acesPoints1;
    std::chrono::system_clock::time_point timeStart, timeStop;
    double sigma = 0.05;
    int fourierOrder = 20;

    ars1.setARSFOrder(fourierOrder);
    ars2.setARSFOrder(fourierOrder);

    rangeToPoint(acesRanges1, 180, -0.5 * M_PI, M_PI / 180.0 * 1.0, acesPoints1);
    acesPoints1.push_back(ars::Point2::Zero());
    std::cout << "Number of input points: " << acesPoints1.size() << std::endl;
    //    for (int i = 0; i < acesPoints1.size(); ++i) {
    //        std::cout << i << "\t" << acesPoints1[i].x() << "\t" << acesPoints1[i].y() << std::endl;
    //    }

    const int trialNum = 1;
    ars1.initLUT(0.0001);
    ars1.setComputeMode(ars::AngularRadonSpectrum2d::PNEBI_DOWNWARD);
    for (int threadNum = 1; threadNum <= 1; ++threadNum) {
        ars1.setThreadNumOMP(threadNum);
        timeStart = std::chrono::system_clock::now();
        for (int t = 0; t < trialNum; ++t) {
            ars1.insertIsotropicGaussians(acesPoints1, sigma);

            std::cout << "trial " << t << ": ars.coefficients().at(0) " << ars1.coefficients().at(0) << ", ars.coefficients().at(2) " << ars1.coefficients().at(2) << std::endl;
        }
        timeStop = std::chrono::system_clock::now();
        double timeAvg = (double) std::chrono::duration_cast<std::chrono::milliseconds>(timeStop - timeStart).count() / (double) trialNum;
        std::cout << "insertIsotropicGaussians() [" << threadNum << " threads]: " << timeAvg << " ms" << std::endl;
    }

    std::cout << "\n------\n" << std::endl;

    ars2.setComputeMode(ars::AngularRadonSpectrum2d::PNEBI_LUT);
    for (int threadNum = 1; threadNum <= 1; ++threadNum) {
        ars2.setThreadNumOMP(threadNum);
        timeStart = std::chrono::system_clock::now();
        for (int t = 0; t < trialNum; ++t) {
            ars2.insertIsotropicGaussians(acesPoints1, sigma);
            std::cout << "trial " << t << ": ars2.coefficients().at(0) " << ars2.coefficients().at(0) << ", ars2.coefficients().at(2) " << ars2.coefficients().at(2) << std::endl;
        }
        timeStop = std::chrono::system_clock::now();
        double timeAvg = (double) std::chrono::duration_cast<std::chrono::milliseconds>(timeStop - timeStart).count() / trialNum;
        std::cout << "insertIsotropicGaussians() [" << threadNum << " threads]: " << timeAvg << " ms" << std::endl;
    }

    std::cout << "\nARS Coefficients:\n";
    std::cout << "\ti \tDownward \tLUT\n";
    for (int i = 0; i < ars1.coefficients().size() && i < ars2.coefficients().size(); ++i) {
        std::cout << "\t" << i << " \t" << ars1.coefficients().at(i) << " \t" << ars2.coefficients().at(i) << "\n";
    }
    std::cout << std::endl;

    std::vector<double> funcFourierRecursDownLUT;
    std::vector<double> funcFourierRecursDown;
    int thnum = 360;
    double dtheta = M_PI / thnum;
    double theta;
    for (int i = 0; i < thnum; ++i) {
        theta = dtheta * i;
        funcFourierRecursDownLUT.push_back(ars1.eval(theta));
        funcFourierRecursDown.push_back(ars2.eval(theta));
    }

    std::cout << "\nBranch and Bound limits:\n";
    int bbnum = 32;
    std::vector<BoundInterval> bbbs(bbnum);
    for (int i = 0; i < bbnum; ++i) {
        bbbs[i].x0 = M_PI * i / bbnum;
        bbbs[i].x1 = M_PI * (i + 1) / bbnum;
        //emotion::findARSFLU(ars.coeffsRecursDown(),bbbs[i].x0,bbbs[i].x1,bbbs[i].y0,bbbs[i].y1);
        ars::findLUFourier(ars1.coefficients(), bbbs[i].x0, bbbs[i].x1, bbbs[i].y0, bbbs[i].y1);
        std::cout << i << ": x0 " << RAD2DEG(bbbs[i].x0) << " x1 " << RAD2DEG(bbbs[i].x1) << ", y0 " << bbbs[i].y0 << " y1 " << bbbs[i].y1 << std::endl;
    }
    
    ars::FourierOptimizerBB1D optim(ars1.coefficients());
    double xopt, ymin, ymax;
    optim.enableXTolerance(true);
    optim.enableYTolerance(true);
    optim.setXTolerance(M_PI/180.0*0.5);
    optim.setYTolerance(1.0);
    optim.findGlobalMax(0,M_PI,xopt,ymin,ymax);
    std::cout << "\n****\nMaximum in x = " << xopt << " (" << RAD2DEG(xopt) << " deg), maximum between [" << ymin << "," << ymax << "]" << std::endl;


    Gnuplot gp("gnuplot -persist");
    double vieweps = 5e-3;
    //  std::ostream& gp = std::cout;
    gp << "set term wxt 0\n";
    gp << "plot '-' title \"ars\" w l, '-' title \"ars lut\" w l, '-' title \"bb\" w l\n";
    for (int i = 0; i < thnum; ++i) {
        gp << (180.0 / thnum * i) << " " << (funcFourierRecursDown[i]) << "\n";
    }
    gp << "e" << std::endl;
    for (int i = 0; i < thnum; ++i) {
        gp << (180.0 / thnum * i) << " " << (funcFourierRecursDownLUT[i]) << "\n";
    }
    gp << "e" << std::endl;
    plotBranchBoundBox(gp, bbbs);
    gp << "e" << std::endl;
    //  std::cout << "\nOptimization: " << std::endl;
    //  double thetaMax, arsfMax;
    //  emotion::findARSFMaxBB(ars.coeffsRecursDown(),0,M_PI,M_PI/180.0*0.05,10.0,thetaMax,arsfMax);
    //  std::cout << "\nMaximum in " << (180.0/M_PI*thetaMax) << " deg, max value " << arsfMax << std::endl;

    return 0;
}

void rangeToPoint(double* ranges, int num, double angleMin, double angleRes, ars::Point2Vector& points) {
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


