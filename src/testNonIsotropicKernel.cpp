#include <iostream>
#include <ars/NonIsotropicKernel.h>
#include <ars/definitions.h>
#include <ars/functions.h>
#include <ars/thirdparty/gnuplot-iostream.h>
#include <ars/Profiler.h>

const double PLOT_EPS = 1e-4;

int main(int argc, char** argv) {
    ars::NonIsotropicKernel nik;
    ars::Vector2 mean1, mean2, mean12, v;
    ars::Matrix2 covar1, covar2, covar12;
    std::vector<double> kernelValRaw;
    std::vector<double> kernelValPolar;
    std::vector<double> kernelValFourier;
    std::vector<double> varianceValRaw;
    std::vector<double> varianceValCos;
    std::vector<double> meanSqValRaw;
    std::vector<double> meanSqValCos;
    std::vector<double> fourierCoeffs;
    double t, fmuDirect, fvarDirect;
    int tnum = 360, fnum = 20;

    // Vectors containing the sampled values of different evaluations of kernels
    kernelValRaw.resize(tnum + 1);
    kernelValPolar.resize(tnum + 1);
    kernelValFourier.resize(tnum + 1);
    meanSqValRaw.resize(tnum + 1);
    meanSqValCos.resize(tnum + 1);
    varianceValRaw.resize(tnum + 1);
    varianceValCos.resize(tnum + 1);


    mean1 << 1.0, 0.5;
    mean2 << -1.0, 3.0;
    covar1 << 1.7, -0.4,
            -0.4, 0.8;
    covar2 << 0.9, 0.2,
            0.2, 2.1;
    mean12 = mean2 - mean1;
    covar12 = covar1 + covar2;

    std::cout << "Gaussian1: mu " << mean1.transpose() << "\ncovar:\n" << covar1 << "\n"
            << "Gaussian2: mu " << mean2.transpose() << "\ncovar:\n" << covar2 << "\n"
            << "Sum covariance matrix:\n" << covar12 << "\n";

    nik.init(mean1, covar1, mean2, covar2);
    nik.computeFourier(fnum, 720, fourierCoeffs);
    
    std::cout << "Fourier coefficients:\n";
    for (int f = 0; f < fnum; ++f) {
        std::cout << "  " << f << " cos: " << fourierCoeffs[2*f] << "\n"
                << "  " << f << " sin: " << fourierCoeffs[2*f+1] << "\n";
    }
    std::cout << std::endl;


    // Computes the raw evalutaion and in polar form
    for (int it = 0; it <= tnum; ++it) {
        t = 2.0 * M_PI * it / tnum;

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
    gp << "set term wxt 1\n";
    gp << "plot '-' title \"raw\" w l, '-' title \"polar\" w l, '-' title \"fourier\" w l"
//       << "'-' title \"variance raw\" w l, '-' title \"variance cos\" w l, "
//       << "'-' title \"mean raw\" w l, '-' title \"mean cos\" w l"
       << "\n";
    for (int it = 0; it <= tnum; ++it) {
        t = 360.0 * it / tnum;
        gp << t << " " << kernelValRaw[it] << "\n";
    }
    gp << "e\n";
    for (int it = 0; it <= tnum; ++it) {
        t = 360.0 * it / tnum;
        gp << t << " " << kernelValPolar[it] << "\n";
    }
    gp << "e\n";
    for (int it = 0; it <= tnum; ++it) {
        t = 360.0 * it / tnum;
        gp << t << " " << kernelValFourier[it] << "\n";
    }
    gp << "e\n";
//    for (int it = 0; it <= tnum; ++it) {
//        t = 360.0 * it / tnum;
//        gp << t << " " << varianceValRaw[it] << "\n";
//    }
//    gp << "e\n";
//    for (int it = 0; it <= tnum; ++it) {
//        t = 360.0 * it / tnum;
//        gp << t << " " << (varianceValCos[it] + PLOT_EPS) << "\n";
//    }
//    gp << "e\n";
//    for (int it = 0; it <= tnum; ++it) {
//        t = 360.0 * it / tnum;
//        gp << t << " " << meanSqValRaw[it] << "\n";
//    }
//    gp << "e\n";
//    for (int it = 0; it <= tnum; ++it) {
//        t = 360.0 * it / tnum;
//        gp << t << " " << (meanSqValCos[it] + PLOT_EPS) << "\n";
//    }
//    gp << "e\n";


    return 0;
}

