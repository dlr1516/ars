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

#include <eigen3/Eigen/Dense>

#include <ars/definitions.h>
#include <ars/ars2d.h>
#include <ars/BBOptimizer1d.h>
#include <ars/HoughSpectrum.h>
#include <ars/HistogramCircularCorrelation.h>

#include <ars/ParamMap.h>
#include <ars/GaussianMixtureEstimator.h>
#include "ars/Profiler.h"
#include "ars/utils.h"



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


int readPoints(std::string filename, ars::VectorVector2 &points);

void rangeToPoint(double* ranges, int num, double angleMin, double angleRes, ars::VectorVector2& points);

void plotBranchBoundBox(std::ostream& out, const std::vector<BoundInterval>& bbbs);

int main(int argc, char **argv) {
    std::string filenameSrc, filenameDst, filenameCfg;
    ars::AngularRadonSpectrum2d arsSrc;
    ars::AngularRadonSpectrum2d arsDst;
    ars::HoughSpectrum hs1;
    ars::HoughSpectrum hs2;
    ars::HistogramCircularCorrelation corr;
    ars::VectorVector2 pointsSrc;
    ars::VectorVector2 pointsDst;
    std::vector<double> correlationFourier;
    std::chrono::system_clock::time_point timeStart, timeStop;
    double sigma = 0.05;
    int fourierOrder = 20;
    double thetaMax1, arsMax1, thetaMax2, arsMax2;
    double thetaTol, fourierTol, thetaMax, corrMax;
    int thnum = 360;
    int hsShiftMax;
    double hsCorrMax;

    std::string clusterAlg;
    double distanceGap, distanceSplit, sigmaMin, weightSum, lmin, lmax, theta, th;
    double clusterDist, meanShiftTol, chi2conf, iseThresh, gaussRes, covarWidth, inlierPerc;
    int arsOrder;
    double arsSigma, arsThetaToll;
    double rotTrue, rotArs;


    ars::ParamMap params;


    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, "");
    std::cout << "config filename: " << filenameCfg << std::endl;
    if (filenameCfg != "") {
        params.read(filenameCfg);
    }

    params.getParam<std::string>("src", filenameSrc, boost::filesystem::current_path().string() + "/points1.txt");
    params.getParam<std::string>("dst", filenameDst, boost::filesystem::current_path().string() + "/points2.txt");
    params.getParam<double>("distanceGap", distanceGap, double(0.6));
    params.getParam<double>("distanceSplit", distanceSplit, double(0.2));
    params.getParam<double>("sigmaMin", sigmaMin, double(0.5));
    params.getParam<double>("clusterDist", clusterDist, double(4.0));
    params.getParam<double>("meanShiftTol", meanShiftTol, double(2.0));
    params.getParam<double>("chi2conf", chi2conf, double(0.80));
    params.getParam<double>("iseThresh", iseThresh, double(0.03));
    params.getParam<double>("gaussRes", gaussRes, double(1.0));
    params.getParam<double>("covarWidth", covarWidth, double(0.2));
    params.getParam<double>("inlierPerc", inlierPerc, double(0.60));
    params.getParam<int>("arsOrder", arsOrder, 20);
    params.getParam<double>("arsSigma", arsSigma, 1.0);
    params.getParam<double>("arsTollDeg", arsThetaToll, 1.0);
    arsThetaToll *= M_PI / 180.0;
    params.getParam<std::string>("clusterAlg", clusterAlg, std::string("hier"));


    ars::GaussianMixtureEstimatorScan *gmeScan = nullptr;
    ars::GaussianMixtureEstimatorHierarchical* gmeHier = nullptr;
    ars::GaussianMixtureEstimatorMeanShift* gmeMean = nullptr;
    ars::GaussianMixtureEstimator* gme = nullptr;

    arsSrc.setARSFOrder(fourierOrder);
    arsDst.setARSFOrder(fourierOrder);

    // Tries to read points from file: if not, the points are read from the available example
    if (readPoints(filenameSrc, pointsSrc) > 0 && readPoints(filenameDst, pointsDst) > 0) {
        std::cout << "Read points from files \"" << filenameSrc << "\" and \"" << filenameDst << "\""
                << std::endl;
    } else {
        std::cout << "Default scan points" << std::endl;


        rangeToPoint(acesRanges1, 180, -0.5 * M_PI, M_PI / 180.0 * 1.0, pointsSrc);
        pointsSrc.push_back(ars::Vector2::Zero());
        std::cout << "Number of input points in scan 1: " << pointsSrc.size() << std::endl;

        rangeToPoint(acesRanges2, 180, -0.5 * M_PI, M_PI / 180.0 * 1.0, pointsDst);
        pointsDst.push_back(ars::Vector2::Zero());
        std::cout << "Number of input points in scan 2: " << pointsDst.size() << std::endl;
    }


    std::cout << "Source size " << pointsSrc.size() << " ; Dest size " << pointsDst.size() << std::endl;


    if (clusterAlg == "scan") {
        gmeScan = new ars::GaussianMixtureEstimatorScan;
        gmeScan->setDistanceGap(distanceGap);
        gmeScan->setDistanceSplit(distanceSplit);
        gmeScan->setSigmaMin(sigmaMin);
        gme = gmeScan;
    } else if (clusterAlg == "hier") {
        gmeHier = new ars::GaussianMixtureEstimatorHierarchical;
        gmeHier->setCellSizeMax(gaussRes);
        gmeHier->setSigmaMin(sigmaMin);
        gmeHier->setIseThreshold(iseThresh);
        gme = gmeHier;
    } else {
        gmeMean = new ars::GaussianMixtureEstimatorMeanShift;
        gmeMean->setSigmaMin(sigmaMin);
        gme = gmeMean;
    }

    // Computes NiArs on src
    std::vector<double> coeffsSrc, coeffsDst, coeffsCor;
    //    double thetaMax, corrMax, fourierTol;

    fourierTol = 1.0; // TODO: check for a proper tolerance

    ars::VectorVector2 means;
    ars::VectorMatrix2 covars;
    std::vector<double> weights;



    // Computes NiArs on src

    gme->clearGaussians();
    {
        ars::ScopedTimer timerGme("GaussianMixtureEstimator::compute()");
        gme->compute(pointsSrc);
    }
    std::cout << "\nFound GMM with " << gme->size() << " kernels:\n";
    weightSum = 0.0;
    for (int i = 0; i < gme->size(); ++i) {
        //            ars::diagonalize(gme->covariance(i), lmin, lmax, theta);
        //            std::cout << "---\n " << i << ": weight " << gme->weight(i) << ", "
        //                    << "mean [" << gme->mean(i).transpose() << "], covar\n"
        //                    << gme->covariance(i) << "\n"
        //                    << "  (lmin " << lmin << ", lmax " << lmax << ", theta[deg] " << (180.0 / M_PI * theta) << ")\n";
        weightSum += gme->weight(i);
    }
    std::cout << "***\nweight sum: " << weightSum << std::endl;
    gme->exportGaussians(means, covars, weights);


    //    {
    //        ars::ScopedTimer("AngularRadonSpectrum2d::insertAnisotropicGaussians()");
    //        if (keepIso) {
    //            setupAnisoAsIso(pointsSrc, means, covars, weights, arsSigma);
    //            arsSrc.insertAnisotropicGaussians(means, covars, weights);
    //        } else
    //            arsSrc.insertAnisotropicGaussians(means, covars, weights);
    //    }

    std::cout << "insertAnisotropicGaussians() src data: means.size() " << means.size() << std::endl;
    arsSrc.setARSFOrder(arsOrder);
    //arsSrc.setAnisotropicStep(arsStep);
    {
        ars::ScopedTimer timer("AngularRadonSpectrum2d::insertAnisotropicGaussians()");
        arsSrc.insertAnisotropicGaussians(means, covars, weights);
        std::cout << "insertAnisotropicGaussians() timer: " << timer.elapsedTimeMs() << " ms" << std::endl;
    }
    coeffsSrc = arsSrc.coefficients();



    //Compute NiArs on dst

    gme->clearGaussians();

    {
        ars::ScopedTimer timerGme("GaussianMixtureEstimator::compute()");
        gme->compute(pointsDst);
    }

    std::cout << "\nFound GMM with " << gme->size() << " kernels:\n";
    weightSum = 0.0;
    for (int i = 0; i < gme->size(); ++i) {
        ars::diagonalize(gme->covariance(i), lmin, lmax, theta);
        //            std::cout << "---\n " << i << ": weight " << gme->weight(i) << ", "
        //                    << "mean [" << gme->mean(i).transpose() << "], covar\n"
        //                    << gme->covariance(i) << "\n"
        //                    << "  (lmin " << lmin << ", lmax " << lmax << ", theta[deg] " << (180.0 / M_PI * theta) << ")\n";

        weightSum += gme->weight(i);
    }
    std::cout << "***\nweight sum: " << weightSum << std::endl;
    gme->exportGaussians(means, covars, weights);


    std::cout << "insertAnisotropicGaussians() dst data: means.size() " << means.size() << std::endl;
    arsDst.setARSFOrder(arsOrder);
    {
        ars::ScopedTimer timer("AngularRadonSpectrum2d::insertAnisotropicGaussians()");
        arsDst.insertAnisotropicGaussians(means, covars, weights);
        std::cout << "insertAnisotropicGaussians() timer: " << timer.elapsedTimeMs() << " ms" << std::endl;
    }
    coeffsDst = arsDst.coefficients();

    //    {
    //        ars::ScopedTimer("AngularRadonSpectrum2d::insertAnisotropicGaussians()");
    //
    //        if (keepIso) {
    //            setupAnisoAsIso(pointsDst, means, covars, weights, arsSigma);
    //            arsDst.insertAnisotropicGaussians(means, covars, weights);
    //        } else
    //            arsDst.insertAnisotropicGaussians(means, covars, weights);
    //    }



    //compute correlation between the two
    {
        ars::ScopedTimer("ars correlation");

        ars::computeFourierCorr(arsSrc.coefficients(), arsDst.coefficients(), correlationFourier);


        ars::findGlobalMaxBBFourier(correlationFourier, 0.0, M_PI, arsThetaToll, fourierTol, thetaMax, corrMax);
        rotArs = thetaMax;
    }

    std::cout << "ARS: best correlation for rotation " << (180.0 / M_PI * thetaMax) << " [deg] with max value " << corrMax << std::endl;

    //    hs1.init(M_PI / thnum, 0.20, 30.0);
    //    hs2.init(M_PI / thnum, 0.20, 30.0);
    //    hs1.insertPoint(pointsSrc.begin(), pointsSrc.end());
    //    hs2.insertPoint(pointsDst.begin(), pointsDst.end());
    //    corr.computeHistSimpleShift(hs1.spectrum(), hs2.spectrum(), thnum, hsShiftMax, hsCorrMax);
    //    std::cout << "HS: best correlation for rotation " << (180.0 * hsShiftMax / thnum) << " [deg] with max value " << corrMax << std::endl;


    Gnuplot gp("gnuplot -persist");
    //    double vieweps = 5e-3;
    //    //  std::ostream& gp = std::cout;
    gp << "set term wxt 0\n";
    //gp << "set terminal postscript eps enhanced color \"Times-Roman\" 24\n";
    //gp << "set output \"registration_ars.eps\"\n";
    gp << "plot '-' title \"ARS 1\" w l lw 5.0, '-' title \"ARS 2\" w l lw 2.0, '-' title \"ARS corr\" w l lw 5.0\n";
    arsMax1 = arsSrc.findMax(thetaMax1);
    for (int i = 0; i < thnum; ++i) {
        double a = M_PI / thnum * i;
        gp << (180.0 / thnum * i) << " " << arsSrc.eval(a) / arsMax1 << "\n";
    }
    gp << "e" << std::endl;
    arsMax2 = arsDst.findMax(thetaMax2);
    for (int i = 0; i < thnum; ++i) {
        double a = M_PI / thnum * i;
        gp << (180.0 / thnum * i) << " " << arsDst.eval(a) / arsMax2 << "\n";
    }
    gp << "e" << std::endl;
    for (int i = 0; i < thnum; ++i) {
        double a = M_PI / thnum * i;
        gp << (180.0 / thnum * i) << " " << ars::evaluateFourier(correlationFourier, 2.0 * a) / corrMax << "\n";
    }
    gp << "e" << std::endl;

#ifdef PLOT_ALL_DEFINED
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
    for (auto& p : pointsSrc) {
        gp << p.x() << " " << p.y() << "\n";
    }
    gp << "e" << std::endl;
    for (auto& p : pointsDst) {
        gp << p.x() << " " << p.y() << "\n";
    }
    gp << "e" << std::endl;
#endif

    return 0;
}

int readPoints(std::string filename, ars::VectorVector2 &points) {
    std::string line, comment;
    ars::Vector2 p;
    size_t pos;
    int count;

    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Cannot open file \"" << filename << "\"" << std::endl;
        return 0;
    }

    points.clear();
    count = 0;
    while (!file.eof()) {
        std::getline(file, line);
        // Remove comments starting with '#'
        comment = "";
        pos = line.find_first_of('#');
        if (pos != std::string::npos) {
            comment = line.substr(pos + 1, line.size());
            line = line.substr(0, pos);
        }
        // Parse the line (after comment removal
        std::stringstream ssline(line);
        if (ssline >> p.x() >> p.y()) {
            points.push_back(p);
            count++;
        }
    }
    file.close();

    return count;
}

void rangeToPoint(double* ranges, int num, double angleMin, double angleRes, ars::VectorVector2& points) {
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


