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
#include <fstream>
#include <Eigen/Dense>

#include <ars/definitions.h>
#include <ars/ars2d.h>
#include <ars/BBOptimizer1d.h>
#include <ars/HoughSpectrum.h>
#include <ars/HistogramCircularCorrelation.h>

#include <ars/ParamMap.h>
#include <ars/GaussianMixtureEstimator.h>
#include <ars/Profiler.h>
#include <ars/utils.h>
#include <ars/HoughTranslationEstimator.h>

#include <chrono>
#include <ars/thirdparty/gnuplot-iostream.h>

#define PRINT_DIM(X) std::cout << #X << " rows " << X.rows() << " cols " << X.cols() << std::endl;
#define RAD2DEG(X) (180.0/M_PI*(X))

int csvToArsVector(const std::string&, ars::VectorVector2&);
int csvToCoeffs(const std::string&, std::vector<double>&);

int main(int argc, char **argv) {
    std::string filenameSrc, filenameDst;
    std::string filenameCoeffSrc, filenameCoeffDst;
    ars::VectorVector2 pointsSrc;
    ars::VectorVector2 pointsDst;
    std::vector<double> coeffSrc;
    std::vector<double> coeffDst;
    std::vector<double> correlationFourier;
    std::chrono::system_clock::time_point timeStart, timeStop;
    double fourierTol, thetaMax, corrMax;
    int thnum = 360;

    std::string clusterAlg;
    double arsThetaToll;
    double rotArs, rotArs2;

    fourierTol = 1.0;
    arsThetaToll = M_PI / 180.0;

    if (argc < 4){
        std::cout << "Use: ars_tester <src cloud> <src coeffs> <dst cloud> <dst coeffs>" <<
            std::endl;
        return -1;
    }
    filenameSrc = argv[1];
    filenameDst = argv[3];
    filenameCoeffSrc = argv[2];
    filenameCoeffDst = argv[4];

    // Tries to read points from files
    if (csvToArsVector(filenameSrc, pointsSrc) != -1) {
        std::cout << "Read points from file " << filenameSrc << std::endl;
    } else {
        std::cout << "Invalid source cloud file!" << std::endl;
        return -1;
    }
    if (csvToArsVector(filenameDst, pointsDst) != -1) {
        std::cout << "Read points from file " << filenameDst << std::endl;
    } else {
        std::cout << "Invalid destination cloud file!" << std::endl;
        return -1;
    }

    std::cout << "Source size " << pointsSrc.size() << " ; Dest size " << pointsDst.size() << std::endl;

    // Tries to read coefficients from files
    if (csvToCoeffs(filenameCoeffSrc, coeffSrc) != -1) {
        std::cout << "Read coefficients from file " << filenameCoeffSrc << std::endl;
    } else {
        std::cout << "Invalid source coefficients file!" << std::endl;
        return -1;
    }
    if (csvToCoeffs(filenameCoeffDst, coeffDst) != -1) {
        std::cout << "Read coefficients from file " << filenameCoeffDst << std::endl;
    } else {
        std::cout << "Invalid destination coefficients file!" << std::endl;
        return -1;
    }
    if(coeffSrc.size() != coeffDst.size()){
        std::cout << "Fourier coefficients have different depths!" << std::endl;
        return -1;
    }
    else
        std::cout << "Fourier coefficients depth: " << (coeffSrc.size()/2) -1  << std::endl;

    //compute correlation between the two
    {
        ars::ScopedTimer("ars correlation");

        ars::computeFourierCorr(coeffSrc, coeffDst, correlationFourier);

        ars::findGlobalMaxBBFourier(correlationFourier, 0.0, M_PI, arsThetaToll, fourierTol, thetaMax, corrMax);
        rotArs = thetaMax;
        rotArs2 = rotArs + M_PI;
    }

    std::cout << "ARS: best correlation for rotation " << (180.0 / M_PI * thetaMax) << 
        " [deg] with max value " << corrMax << std::endl;

    //plot dst, src and transformed src
    ars::VectorVector2 pointsTrans;
    for(auto& p : pointsSrc){
        ars::Vector2 tmp;
        tmp.x() = (p.x()*std::cos(rotArs))-(p.y()*std::sin(rotArs));
        tmp.y() = (p.x()*std::sin(rotArs))+(p.y()*std::cos(rotArs));
        pointsTrans.push_back(tmp);
    }

    ars::VectorVector2 pointsTrans2;
    for(auto& p : pointsSrc){
        ars::Vector2 tmp;
        tmp.x() = (p.x()*std::cos(rotArs2))-(p.y()*std::sin(rotArs2));
        tmp.y() = (p.x()*std::sin(rotArs2))+(p.y()*std::cos(rotArs2));
        pointsTrans2.push_back(tmp);
    }

    timeStart = std::chrono::system_clock::now();
    ars::HoughTranslationEstimator::Translation tx, ty, tx2, ty2;
    ars::HoughTranslationEstimator tEst;
    double res = 0.01;
    double range = 20.0;
    double sigma = 0.05;
    tEst.init(res, range, sigma);

    tEst.compute(pointsTrans, pointsDst, tx, ty);
    timeStop = std::chrono::system_clock::now();
    double time = (double) std::chrono::duration_cast<std::chrono::milliseconds>(timeStop - timeStart).count();
    {
        Gnuplot gp("gnuplot -persist");
        gp << "set term wxt 0 title 'X histograms 1'\n";
        gp << "set style fill transparent solid 0.7\n";
        gp << "plot '-' title 'src' with boxes,'-' title 'dst' with boxes\n";
        tEst.plotXHistograms(gp);
    }
    {
        Gnuplot gp("gnuplot -persist");
        gp << "set term wxt 0  title 'Y histograms 1'\n";
        gp << "set style fill transparent solid 0.7\n";
        gp << "plot '-' title 'src' with boxes,'-' title 'dst' with boxes\n";
        tEst.plotYHistograms(gp);
    }

    timeStart = std::chrono::system_clock::now();
    tEst.compute(pointsTrans2, pointsDst, tx2, ty2);
    timeStop = std::chrono::system_clock::now();
    time += (double) std::chrono::duration_cast<std::chrono::milliseconds>(timeStop - timeStart).count();
    std::cout << "translation time: " << time << " ms" << std::endl;
    {
        Gnuplot gp("gnuplot -persist");
        gp << "set term wxt 0 title 'X histograms 2'\n";
        gp << "set style fill transparent solid 0.7\n";
        gp << "plot '-' title 'src' with boxes,'-' title 'dst' with boxes\n";
        tEst.plotXHistograms(gp);
    }
    {
        Gnuplot gp("gnuplot -persist");
        gp << "set term wxt 0  title 'Y histograms 2'\n";
        gp << "set style fill transparent solid 0.7\n";
        gp << "plot '-' title 'src' with boxes,'-' title 'dst' with boxes\n";
        tEst.plotYHistograms(gp);
    }

    std::cout << "translation 1: " << tx.translation << ", " << ty.translation << 
        " with confidence " << tx.correlation + ty.correlation << std::endl;
    std::cout << "translation 2: " << tx2.translation << ", " << ty2.translation << 
        " with confidence " << tx2.correlation + ty2.correlation << std::endl;

    ars::VectorVector2 pointsTransl;
    for(auto& p : pointsTrans){
        ars::Vector2 tmp;
        tmp.x() = p.x() - tx.translation;
        tmp.y() = p.y() - ty.translation;
        pointsTransl.push_back(tmp);
    }

    ars::VectorVector2 pointsTransl2;
    for(auto& p : pointsTrans2){
        ars::Vector2 tmp;
        tmp.x() = p.x() - tx2.translation;
        tmp.y() = p.y() - ty2.translation;
        pointsTransl2.push_back(tmp);
    }

    {
        Gnuplot gp("gnuplot -persist");
        gp << "set term wxt 0\n";
        gp << "set size ratio -1\n";
        gp << "plot '-' title \"src\" w p pt 7 ps 0.7, '-' title \"dst\" w p pt 7 ps 0.7," << 
            "'-' title \"transformed\" w p pt 7 ps 0.5, '-' title \"transformed 2\" w p pt 7 ps 0.5\n";

        for (auto& point : pointsSrc) {
            gp << point.x() << " " << point.y() << "\n";
        }
        gp << "e" << std::endl;

        for (auto& point : pointsDst) {
            gp << point.x() << " " << point.y() << "\n";
        }
        gp << "e" << std::endl;

        for (auto& point : pointsTrans) {
            gp << point.x() << " " << point.y() << "\n";
        }
        gp << "e" << std::endl;

        for (auto& point : pointsTrans2) {
            gp << point.x() << " " << point.y() << "\n";
        }
        gp << "e" << std::endl;
    }
    {
        Gnuplot gp("gnuplot -persist");
        gp << "set term wxt 0\n";
        gp << "set size ratio -1\n";
        gp << "plot '-' title \"dst\" w p pt 7 ps 0.7, '-' title \"translated\" w p pt 7 ps 0.7," << 
            "'-' title \"translated 2\" w p pt 7 ps 0.5\n";

        for (auto& point : pointsDst) {
            gp << point.x() << " " << point.y() << "\n";
        }
        gp << "e" << std::endl;

        for (auto& point : pointsTransl) {
            gp << point.x() << " " << point.y() << "\n";
        }
        gp << "e" << std::endl;

        for (auto& point : pointsTransl2) {
            gp << point.x() << " " << point.y() << "\n";
        }
        gp << "e" << std::endl;
    }

    return 0;
}

int csvToArsVector(const std::string& filename, ars::VectorVector2& cloud){
    std::ifstream file(filename);

    std::string line, comment;
    ars::Vector2 p;
    size_t pos;

    if (!file) {
        std::cerr << "Cannot open file \"" << filename << "\"" << std::endl;
        return -1;
    }
    cloud.clear();
    std::getline(file, line); //discard first line of csv
    while (!file.eof()) {
        std::getline(file, line);
        // Remove comments starting with '#'
        comment = "";
        pos = line.find_first_of('#');
        if (pos != std::string::npos) {
            comment = line.substr(pos + 1, line.size());
            line = line.substr(0, pos);
        }
        // Parse the line (after comment removal)
        pos = line.find(',');
        if (pos != std::string::npos) {
            line.replace(pos, 1, " ");
        }
        std::stringstream ssline(line);
        if (ssline >> p.x() >> p.y()) {
            cloud.push_back(p);
        }
    }
    file.close();
    return 0;
}

int csvToCoeffs(const std::string& filename, std::vector<double>& coeffs){
    std::ifstream file(filename);
    std::string line, comment;
    coeffs.clear();
    while (!file.eof()) {
        double tmp = 0;
        std::getline(file, line, ',');
        std::stringstream ssline(line);
        if (ssline >> tmp) {
            coeffs.push_back(tmp);
        }
    }
    file.close();
    if(coeffs.size() > 0)
        return 0;
    else
        return -1;
}