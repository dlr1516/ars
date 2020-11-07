#include <iostream>
#include <ars/AnisotropicKernel.h>
#include <ars/ars2d.h>
#include <ars/definitions.h>
#include <ars/functions.h>
#include <ars/thirdparty/gnuplot-iostream.h>
#include <ars/Profiler.h>
#include <ars/ParamMap.h>

const double PLOT_EPS = 1e-4;

int main(int argc, char** argv) {
    ars::AnisotropicKernel nik;
    ars::Vector2 mean1, mean2, mean12, v;
    ars::Matrix2 covar1, covar2, covar12;
    std::vector<double> kernelValRaw;
    std::vector<double> kernelValPolar;
    std::vector<double> kernelValFourier;
    std::vector<double> varianceValRaw;
    std::vector<double> varianceValCos;
    std::vector<double> meanSqValRaw;
    std::vector<double> meanSqValCos;
    std::vector<double> fourierCoeffs, fourierCoeffsIsot;
    double t, fmuDirect, fvarDirect, sigmaVal, sigmaValSq, lambdaSqNorm, phi, normalizer, plotEps;
    int plotStep, arsOrder, arsStep;
    ars::ParamMap params;
    std::string filenameCfg;

    // Reads params from command line
    params.read(argc, argv);
    params.getParam("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    params.getParam<double>("sigmaVal", sigmaVal, double(0.08));
    params.getParam<int>("arsOrder", arsOrder, int(20));
    params.getParam<int>("arsStep", arsStep, int(720));
    params.getParam<int>("plotStep", plotStep, int(1440));
    params.getParam<double>("plotEps", plotEps, double(1e-3));

    std::cout << "\nParams:" << std::endl;
    params.write(std::cout);
    std::cout << "-------\n" << std::endl;

    // Vectors containing the sampled values of different evaluations of kernels
    kernelValRaw.resize(plotStep + 1);
    kernelValPolar.resize(plotStep + 1);
    kernelValFourier.resize(plotStep + 1);
    meanSqValRaw.resize(plotStep + 1);
    meanSqValCos.resize(plotStep + 1);
    varianceValRaw.resize(plotStep + 1);
    varianceValCos.resize(plotStep + 1);


    mean1 << 1.0, 0.5;
    mean2 << -1.0, 3.0;
    covar1 << 1.7, -0.4,
            -0.4, 0.8;
    covar2 << 0.9, 0.2,
            0.2, 2.1;
    mean12 = mean2 - mean1;
    covar12 = covar1 + covar2;

    std::cout << "***\nTEST 1\n"
            << "Gaussian1: mu " << mean1.transpose() << "\ncovar:\n" << covar1 << "\n"
            << "Gaussian2: mu " << mean2.transpose() << "\ncovar:\n" << covar2 << "\n"
            << "Sum covariance matrix:\n" << covar12 << "\n";

    nik.init(mean1, covar1, mean2, covar2);
    nik.computeFourier(arsOrder, arsStep, fourierCoeffs);

    std::cout << "Fourier coefficients:\n";
    for (int f = 0; f < arsOrder; ++f) {
        std::cout << "  " << f << " cos: " << fourierCoeffs[2 * f] << "\n"
                << "  " << f << " sin: " << fourierCoeffs[2 * f + 1] << "\n";
    }
    std::cout << std::endl;

    // Computes the raw evalutaion and in polar form
    for (int it = 0; it <= plotStep; ++it) {
        t = 2.0 * M_PI * it / plotStep;

        // Computes the raw value 
        v << cos(t), sin(t);
        fmuDirect = mean12.dot(v);
        fvarDirect = v.transpose() * covar12 * v;
        kernelValRaw[it] = exp(-0.5 * fmuDirect * fmuDirect / fvarDirect) / sqrt(2.0 * M_PI * fvarDirect);

        // Computes the polar value
        kernelValPolar[it] = nik.value(t);

        // Computes the Fourier value
        kernelValFourier[it] = ars::evaluateFourier(fourierCoeffs, 2.0 * t);

        // Debug partial
        meanSqValRaw[it] = fmuDirect * fmuDirect;
        meanSqValCos[it] = 0.5 * nik.getMuModule() * nik.getMuModule() * (1.0 + cos(2.0 * t - 2.0 * nik.getMuPhase()));
        varianceValRaw[it] = fvarDirect;
        varianceValCos[it] = nik.getVarianceModule() * (1.0 + nik.getVariancePerc() * cos(2.0 * t - 2.0 * nik.getVariancePhase()));
    }

    // Plot
    Gnuplot gp("gnuplot -persist");
    //std::ostream& gp = std::cerr;
    gp << "set term wxt 0 title \"kernel\"\n";
    gp << "plot '-' title \"raw\" w l, '-' title \"polar\" w l, '-' title \"fourier\" w l"
            //       << "'-' title \"variance raw\" w l, '-' title \"variance cos\" w l, "
            //       << "'-' title \"mean raw\" w l, '-' title \"mean cos\" w l"
            << "\n";
    for (int it = 0; it <= plotStep; ++it) {
        t = 360.0 * it / plotStep;
        gp << t << " " << kernelValRaw[it] << "\n";
    }
    gp << "e\n";
    for (int it = 0; it <= plotStep; ++it) {
        t = 360.0 * it / plotStep;
        gp << t << " " << kernelValPolar[it] << "\n";
    }
    gp << "e\n";
    for (int it = 0; it <= plotStep; ++it) {
        t = 360.0 * it / plotStep;
        gp << t << " " << kernelValFourier[it] << "\n";
    }
    gp << "e\n";

    // ----------------------------------------------------

    mean1 << 1.0, 0.5;
    mean2 << -1.0, 3.0;
    sigmaValSq = sigmaVal * sigmaVal;
    covar1 << sigmaValSq, 0.0,
            0.0, sigmaValSq;
    covar2 = covar1;
    mean12 = mean2 - mean1;
    covar12 = covar1 + covar2;
    lambdaSqNorm = 0.125 * mean12.dot(mean12) / sigmaValSq;
    phi = atan2(mean12(1), mean12(0));
    normalizer = 1.0 / sqrt(4.0 * M_PI * sigmaValSq);

    std::cout << "***\nTEST 2: consistency check on isotropic GMM\n"
            << "Gaussian1: mu " << mean1.transpose() << "\ncovar:\n" << covar1 << "\n"
            << "Gaussian2: mu " << mean2.transpose() << "\ncovar:\n" << covar2 << "\n";

    nik.init(mean1, covar1, mean2, covar2);
    nik.computeFourier(arsOrder, arsStep, fourierCoeffs);
    ARS_PRINT("Anisotropic kernel computed parameters:\n"
            << "  muMod_ " << nik.getMuModule()
            << ", muAng_[rad] " << nik.getMuPhase() << " [deg] " << (180.0 / M_PI * nik.getMuPhase()) << "\n"
            << "  sigmaMod_ " << nik.getVarianceModule()
            << ", sigmaAng_[rad] " << nik.getVariancePhase() << " [deg] " << (180.0 / M_PI * nik.getVariancePhase())
            << ", sigmaDif_ " << nik.getVariancePerc()
            << ", sqrt(2 * M_PI * sigmaMod_) " << sqrt(2 * M_PI * nik.getVarianceModule())
            << "\n");

    fourierCoeffsIsot.resize(2 * arsOrder + 2, 0.0);
    ars::updateARSF2CoeffRecursDown(lambdaSqNorm, phi, 1.0, arsOrder, fourierCoeffsIsot);
    ARS_PRINT("Isotropic kernel computed parameters:\n"
            << " 2 * sigmaVal^2 " << (2.0 * sigmaVal * sigmaVal)
            << ", lambda " << sqrt(mean12.dot(mean12))
            << ", lambda^2 / (8 * sigmaVal^2) " << lambdaSqNorm
            << ", phi[rad] " << phi << " [deg] " << (180.0 / M_PI * phi));

    std::cout << "Fourier coefficients (anisotropic - isotropic): "
            << "fourierCoeffs.size() " << fourierCoeffs.size()
            << ", fourierCoeffsIsot.size() " << fourierCoeffsIsot.size() << "\n";
    for (int f = 0; f < arsOrder; ++f) {
        std::cout << "  " << f << " cos: \t" << fourierCoeffs[2 * f] << "\t" << normalizer * fourierCoeffsIsot[2 * f] << "\n"
                << "  " << f << " sin: \t" << fourierCoeffs[2 * f + 1] << "\t" << normalizer * fourierCoeffsIsot[2 * f + 1] << "\n";
    }
    std::cout << std::endl;

    gp << "set term wxt 1 title \"iso-aniso\"\n";
    gp << "plot '-' title \"direct isotr.\" w l, "
            << "'-' title \"direct anisotr.\" w l, "
            << "'-' title \"isotropic\" w l, "
            << "'-' title \"anisotropic\" w l"
            << "\n";
    double ker;
    for (int it = 0; it <= plotStep; ++it) {
        t = M_PI * it / plotStep;
        ker = exp(-lambdaSqNorm * (1.0 + cos(2.0 * t - 2.0 * phi))) / sqrt(4.0 * M_PI * sigmaValSq);
        gp << (180.0 / M_PI * t) << " " << ker << "\n";
    }
    gp << "e\n";
    for (int it = 0; it <= plotStep; ++it) {
        t = M_PI * it / plotStep;
        gp << (180.0 / M_PI * t) << " " << nik.value(t) + plotEps << "\n";
    }
    gp << "e\n";
    for (int it = 0; it <= plotStep; ++it) {
        t = M_PI * it / plotStep;
        gp << (180.0 / M_PI * t) << " " << normalizer * ars::evaluateFourier(fourierCoeffsIsot, 2.0 * t) + 2.0 * plotEps << "\n";
    }
    gp << "e\n";
    for (int it = 0; it <= plotStep; ++it) {
        t = M_PI * it / plotStep;
        gp << (180.0 / M_PI * t) << " " << ars::evaluateFourier(fourierCoeffs, 2.0 * t) + 3.0 * plotEps << "\n";
    }
    gp << "e\n";

    return 0;
}

