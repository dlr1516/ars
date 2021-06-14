#include <iostream>
#include <ars/GaussianMixtureEstimator.h>
#include <ars/ParamMap.h>
#include <ars/thirdparty/gnuplot-iostream.h>
#include <ars/utils.h>

double acesRanges[] = {50.00, 50.00, 50.00, 5.26, 5.21, 5.06, 5.01, 3.01, 2.94, 2.89, 2.84, 2.74, 2.69, 2.64, 2.59, 2.54, 2.49, 2.49, 2.44, 2.39, 2.34, 2.29, 2.29, 2.29, 2.39, 2.39, 2.49, 2.51, 2.61, 2.66, 2.76, 2.81, 2.96, 3.01, 3.11, 3.26, 3.01, 3.01, 3.01, 3.06, 3.21, 6.86, 6.86, 6.81, 6.76, 6.71, 6.71, 6.66, 6.61, 6.66, 6.56, 6.56, 6.56, 6.46, 6.46, 6.41, 6.46, 6.46, 4.11, 3.96, 3.96, 4.96, 4.86, 5.21, 7.41, 4.61, 5.16, 6.26, 6.26, 6.31, 4.86, 5.01, 5.86, 5.81, 4.21, 4.26, 4.31, 4.41, 4.39, 4.46, 5.31, 5.06, 5.26, 4.96, 6.01, 5.76, 5.61, 5.36, 5.26, 5.01, 4.21, 4.16, 4.01, 3.91, 3.61, 3.21, 3.26, 3.16, 3.06, 3.01, 3.31, 3.21, 3.16, 2.16, 2.19, 2.16, 2.21, 2.11, 2.01, 2.01, 2.06, 2.84, 2.91, 2.91, 3.01, 3.11, 3.21, 3.81, 4.06, 7.11, 7.06, 7.01, 6.96, 6.86, 4.31, 6.76, 6.71, 6.66, 6.61, 5.46, 5.41, 6.46, 6.21, 6.31, 6.51, 7.26, 7.46, 50.00, 2.01, 1.94, 1.94, 1.94, 2.31, 1.86, 1.84, 1.84, 1.81, 1.96, 26.46, 20.76, 2.11, 2.12, 2.17, 2.14, 2.09, 2.09, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.14, 2.19, 2.19, 2.24, 2.24, 2.24, 2.24, 2.29, 2.29, 2.29, 2.29, 2.29, 2.39, 2.39, 2.39, 2.44};

void rangeToPoint(double* ranges, int num, double angleMin, double angleRes, ars::VectorVector2& points);

void plotEllipse(std::ostream& out, int idx, const ars::Vector2& mean, const ars::Matrix2& covar);

void plotEllipses(std::ostream& out, const ars::VectorVector2& means, const ars::VectorMatrix2& covars);

int main(int argc, char** argv) {
    ars::VectorVector2 acesPoints, means;
    ars::VectorMatrix2 covars, covarsUniform;
    ars::GaussianMixtureEstimatorScan gmeScan;
    ars::GaussianMixtureEstimatorMeanShift gmeMean;
    double sigmaMin, clusterDist, meanShiftTol, distanceGap, distanceSplit;
    int iterationNumMax;
    ars::ParamMap params;
    std::string filenameCfg;

    // Reads params from command line
    params.read(argc, argv);
    params.getParam("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    params.getParam<double>("sigmaMin", sigmaMin, double(0.05));
    params.getParam<double>("clusterDist", clusterDist, double(4.0));
    params.getParam<double>("meanShiftTol", meanShiftTol, double(2.0));
    params.getParam<int>("iterationNumMax", iterationNumMax, int(30));
    params.getParam<double>("distanceGap", distanceGap, double(0.6));
    params.getParam<double>("distanceSplit", distanceSplit, double(0.2));

    std::cout << "\nParams:" << std::endl;
    params.write(std::cout);
    std::cout << "-------\n" << std::endl;

    rangeToPoint(acesRanges, 180, -0.5 * M_PI, M_PI / 180.0 * 1.0, acesPoints);
    acesPoints.push_back(ars::Vector2::Zero());
    std::cout << "Number of input points: " << acesPoints.size() << std::endl;

    std::cout << "\n---\nTesting GaussianMixtureEstimatorMeanShift:" << std::endl;
    gmeScan.setSigmaMin(sigmaMin);
    gmeScan.setDistanceGap(distanceGap);
    gmeScan.setDistanceSplit(distanceSplit);
    gmeScan.compute(acesPoints);
    std::cout << "Found " << gmeScan.size() << " clusters" << std::endl;

    std::cout << "\n---\nTesting GaussianMixtureEstimatorMeanShift:" << std::endl;
    gmeMean.setSigmaMin(sigmaMin);
    gmeMean.setClusterDistance(clusterDist);
    gmeMean.setMeanShiftTol(meanShiftTol);
    gmeMean.setIterationNumMax(iterationNumMax);
    gmeMean.compute(acesPoints);

    {
        double weightSum, lmin, lmax, theta;
        
        weightSum = 0.0;
        std::cout << "\nFound GMM with " << gmeMean.size() << " kernels:\n";

        for (int i = 0; i < gmeMean.size(); ++i) {
            ars::diagonalize(gmeMean.covariance(i), lmin, lmax, theta);
            std::cout << "---\n " << i << ": weight " << gmeMean.weight(i) << ", "
                    << "mean [" << gmeMean.mean(i).transpose() << "], covar\n"
                    << gmeMean.covariance(i) << "\n"
                    << "  (lmin " << lmin << ", lmax " << lmax << ", theta[deg] " << (180.0 / M_PI * theta) << ")\n";
            weightSum += gmeMean.weight(i);
        }
        std::cout << "***\nweight sum: " << weightSum << std::endl;
    }

    
    Gnuplot gp("gnuplot -persist");
    gp << "set term wxt 0 title \"GMM\"\n";
    gp << "set xrange [0.0:8.0]\n";
    gp << "set yrange [-8.0:8.0]\n";
    gp << "set size ratio -1\n";
    for (int i = 0; i < gmeMean.size(); ++i) {
        plotEllipse(gp, i + 1, gmeMean.mean(i), gmeMean.covariance(i));
    }
    gp << "plot '-' title \"scan\" w p pt 7 ps 0.5\n";
    for (auto& p : acesPoints) {
        gp << p.x() << " " << p.y() << "\n";
    }
    gp << "e\n";

    // Computes ARS

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

