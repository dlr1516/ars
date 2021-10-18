#include <iostream>
#include <ars/GaussianMixtureEstimator.h>
#include <ars/ParamMap.h>
#include <ars/thirdparty/gnuplot-iostream.h>
#include <ars/utils.h>

double acesRanges[] = {50.00, 50.00, 50.00, 5.26, 5.21, 5.06, 5.01, 3.01, 2.94, 2.89, 2.84, 2.74, 2.69, 2.64, 2.59, 2.54, 2.49, 2.49, 2.44, 2.39, 2.34, 2.29, 2.29, 2.29, 2.39, 2.39, 2.49, 2.51, 2.61, 2.66, 2.76, 2.81, 2.96, 3.01, 3.11, 3.26, 3.01, 3.01, 3.01, 3.06, 3.21, 6.86, 6.86, 6.81, 6.76, 6.71, 6.71, 6.66, 6.61, 6.66, 6.56, 6.56, 6.56, 6.46, 6.46, 6.41, 6.46, 6.46, 4.11, 3.96, 3.96, 4.96, 4.86, 5.21, 7.41, 4.61, 5.16, 6.26, 6.26, 6.31, 4.86, 5.01, 5.86, 5.81, 4.21, 4.26, 4.31, 4.41, 4.39, 4.46, 5.31, 5.06, 5.26, 4.96, 6.01, 5.76, 5.61, 5.36, 5.26, 5.01, 4.21, 4.16, 4.01, 3.91, 3.61, 3.21, 3.26, 3.16, 3.06, 3.01, 3.31, 3.21, 3.16, 2.16, 2.19, 2.16, 2.21, 2.11, 2.01, 2.01, 2.06, 2.84, 2.91, 2.91, 3.01, 3.11, 3.21, 3.81, 4.06, 7.11, 7.06, 7.01, 6.96, 6.86, 4.31, 6.76, 6.71, 6.66, 6.61, 5.46, 5.41, 6.46, 6.21, 6.31, 6.51, 7.26, 7.46, 50.00, 2.01, 1.94, 1.94, 1.94, 2.31, 1.86, 1.84, 1.84, 1.81, 1.96, 26.46, 20.76, 2.11, 2.12, 2.17, 2.14, 2.09, 2.09, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.19, 2.19, 2.24, 2.24, 2.24, 2.24, 2.29, 2.29, 2.29, 2.29, 2.29, 2.39, 2.39, 2.39, 2.44};

int readPoints(std::string filename, ars::VectorVector2& points);

void findMinMax(const ars::VectorVector2& points, ars::Vector2& pmin, ars::Vector2& pmax);

void rangeToPoint(double* ranges, int num, double angleMin, double angleRes, double rangeMax, ars::VectorVector2& points);

void plotEllipse(std::ostream& out, int idx, const ars::Vector2& mean, const ars::Matrix2& covar);

void plotEllipses(std::ostream& out, const ars::VectorVector2& means, const ars::VectorMatrix2& covars);

int main(int argc, char** argv) {
    ars::VectorVector2 acesPoints, means;
    ars::VectorMatrix2 covars, covarsUniform;
    ars::GaussianMixtureEstimator* gme = nullptr;
    ars::GaussianMixtureEstimatorScan* gmeScan = nullptr;
    ars::GaussianMixtureEstimatorHierarchical* gmeHier = nullptr;
    ars::GaussianMixtureEstimatorMeanShift* gmeMean = nullptr;
    double sigmaMin, clusterDist, meanShiftTol, distanceGap, distanceSplit, chi2conf, gaussRes, inlierPerc, covarWidth, weightSum, lmin, lmax, theta;
    double iseThresh;
    ars::Vector2 pmin, pmax;
    double rangeMax, framePlot;
    int iterationNumMax;
    ars::ParamMap params;
    std::string filenameCfg, filenameIn, clusterAlg;

    // Reads params from command line
    params.read(argc, argv);
    params.getParam("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    params.getParam<std::string>("in", filenameIn, std::string(""));
    params.getParam<double>("sigmaMin", sigmaMin, double(0.025));
    params.getParam<double>("clusterDist", clusterDist, double(4.0));
    params.getParam<double>("meanShiftTol", meanShiftTol, double(2.0));
    params.getParam<int>("iterationNumMax", iterationNumMax, int(30));
    params.getParam<double>("distanceGap", distanceGap, double(0.6));
    params.getParam<double>("distanceSplit", distanceSplit, double(0.2));
    params.getParam<double>("chi2conf", chi2conf, double(0.80));
    params.getParam<double>("iseThresh", iseThresh, double(0.30));
    params.getParam<double>("gaussRes", gaussRes, double(16.0 * sigmaMin));
    params.getParam<double>("covarWidth", covarWidth, double(0.2));
    params.getParam<double>("inlierPerc", inlierPerc, double(0.60));
    params.getParam<std::string>("clusterAlg", clusterAlg, std::string("scan"));
    params.getParam<double>("rangeMax", rangeMax, double(10.0));

    std::cout << "\nParams:" << std::endl;
    params.write(std::cout);
    std::cout << "-------\n" << std::endl;

    // Tries to read points from file: if not, the points are read from the available example
    if (readPoints(filenameIn, acesPoints) > 0) {
        std::cout << "Read points from file \"" << filenameIn << "\"" << std::endl;
    } else {
        std::cout << "Default scan points" << std::endl;
        rangeToPoint(acesRanges, 180, -0.5 * M_PI, M_PI / 180.0 * 1.0, rangeMax, acesPoints);
        acesPoints.push_back(ars::Vector2::Zero());
    }
    findMinMax(acesPoints, pmin, pmax);
    std::cout << "Number of input points: " << acesPoints.size() << " in region "
            << " pmin [" << pmin.transpose() << "]  pmax [" << pmax.transpose() << "]" << std::endl;
    framePlot = 0.05 * (pmax(0) - pmin(0));
    if (framePlot < 0.05 * (pmax(1) - pmin(1))) {
        framePlot = 0.05 * (pmax(1) - pmin(1));
    }

    if (clusterAlg == "scan") {
        std::cout << "\n---\nTesting GaussianMixtureEstimatorMeanShift:" << std::endl;
        gmeScan = new ars::GaussianMixtureEstimatorScan;
        gmeScan->setDistanceGap(distanceGap);
        gmeScan->setDistanceSplit(distanceSplit);
        gmeScan->setSigmaMin(sigmaMin);
        gme = gmeScan;
    } else if (clusterAlg == "hier") {
        gmeHier = new ars::GaussianMixtureEstimatorHierarchical;
        gmeHier->setSigmaMin(sigmaMin);
        gmeHier->setCovarWidth(covarWidth);
        //        gmeHier->setChiConfidence(chi2conf);
        gmeHier->setIseThreshold(iseThresh);
        gmeHier->setCellSizeMax(gaussRes);
        gme = gmeHier;
    } else {
        std::cout << "\n---\nTesting GaussianMixtureEstimatorMeanShift:" << std::endl;
        gmeMean = new ars::GaussianMixtureEstimatorMeanShift;
        gmeMean->setSigmaMin(sigmaMin);
        gmeMean->setIterationNumMax(iterationNumMax);
        gme = gmeMean;
    }
    gme->compute(acesPoints);
    std::cout << "Found " << gme->size() << " clusters" << std::endl;



    weightSum = 0.0;
    std::cout << "\nFound GMM with " << gme->size() << " kernels:\n";
    for (int i = 0; i < gme->size(); ++i) {
        ars::diagonalize(gme->covariance(i), lmin, lmax, theta);
        //if (gme->size() < 30) {
//        std::cout << "---\n " << i << ": weight " << gme->weight(i) << ", "
//                << "mean [" << gme->mean(i).transpose() << "], covar\n"
//                << gme->covariance(i) << "\n"
//                << "  (lmin " << lmin << ", lmax " << lmax << ", theta[deg] " << (180.0 / M_PI * theta) << ")\n";
        //}
        weightSum += gme->weight(i);
    }
    std::cout << "***\nweight sum: " << weightSum << std::endl;

    Gnuplot gp("gnuplot -persist");
    gp << "set term wxt 0 title \"GMM " << clusterAlg << "\"\n";
    gp << "set xrange [" << (pmin(0) - framePlot) << ":" << (pmax(0) + framePlot) << "]\n";
    gp << "set yrange [" << (pmin(1) - framePlot) << ":" << (pmax(1) + framePlot) << "]\n";
    gp << "set size ratio -1\n";
    for (int i = 0; i < gme->size(); ++i) {
        plotEllipse(gp, i + 1, gme->mean(i), gme->covariance(i));
    }
    gp << "plot '-' title \"scan\" w p pt 7 ps 0.5\n";
    for (auto& p : acesPoints) {
        gp << p.x() << " " << p.y() << "\n";
    }
    gp << "e\n";

    // Computes ARS

    return 0;
}

int readPoints(std::string filename, ars::VectorVector2& points) {
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

void findMinMax(const ars::VectorVector2& points, ars::Vector2& pmin, ars::Vector2& pmax) {
    if (points.empty()) {
        pmin = pmax = ars::Vector2::Zero();
    }
    pmin = points[0];
    pmax = points[0];
    for (int i = 1; i < points.size(); ++i) {
        if (points[i](0) < pmin(0)) pmin(0) = points[i](0);
        if (points[i](0) > pmax(0)) pmax(0) = points[i](0);
        if (points[i](1) < pmin(1)) pmin(1) = points[i](1);
        if (points[i](1) > pmax(1)) pmax(1) = points[i](1);
    }
}

void rangeToPoint(double* ranges, int num, double angleMin, double angleRes, double rangeMax, ars::VectorVector2& points) {
    Eigen::Vector2d p;
    for (int i = 0; i < num; ++i) {
        double a = angleMin + angleRes * i;
        if (rangeMax > 0.0 && ranges[i] < rangeMax) {
            p << ranges[i] * cos(a), ranges[i] * sin(a);
            points.push_back(p);
        }
    }
}

void plotEllipse(std::ostream& out, int idx, const ars::Vector2& mean, const ars::Matrix2& covar) {
    double lmin, lmax, angle;
    //    set object 1 ellipse center 1.5, 1  size 6, 12  angle 60 front fs empty bo 3
    //    plot '-' with points
    // Confidence 0.95 -> chi2 5.991 -> axis 

    ars::diagonalize(covar, lmin, lmax, angle);
    out << "set object " << idx << " ellipse center " << mean(0) << ", " << mean(1)
            << " size " << sqrt(5.991 * lmax) << ", " << sqrt(5.991 * lmin) << " angle " << (180.0 / M_PI * angle)
            << " front fs empty bo 3\n";
}

void plotEllipses(std::ostream& out, const ars::VectorVector2& means, const ars::VectorMatrix2& covars) {
    for (int i = 0; i < means.size() && i < covars.size(); ++i) {
        plotEllipse(out, i, means[i], covars[i]);
    }
}

