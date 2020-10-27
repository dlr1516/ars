#include <iostream>
#include <ars/NonIsotropicKernel.h>
#include <ars/definitions.h>
#include <ars/thirdparty/gnuplot-iostream.h>
#include <ars/Profiler.h>

int main(int argc, char** argv) {
    ars::NonIsotropicKernel nik;
    ars::Vector2 mean1, mean2, mean12, v;
    ars::Matrix2 covar1, covar2, covar12;
    std::vector<double> kernelValRaw;
    std::vector<double> kernelValPolar;
    std::vector<double> kernelValFourier;
    std::vector<double> varianceVal;
    double t, fmuDirect, fvarDirect;
    int tnum = 360;

    // Vectors containing the sampled values of different evaluations of kernels
    kernelValRaw.resize(tnum + 1);
    kernelValPolar.resize(tnum + 1);
    kernelValFourier.resize(tnum + 1);
    varianceVal.push_back(tnum + 1);


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
    
    

    // Computes the raw evalutaion and in polar form
    for (int it = 0; it <= tnum; ++it) {
        t = 2.0 * M_PI * it / tnum;
        
        // Computes the raw value 
        v << cos(t), sin(t);
        fmuDirect = mean12.dot(v);
        fvarDirect = v.transpose() * covar12 * v;
        kernelValRaw[it] = exp(-0.5 * fmuDirect * fmuDirect / fvarDirect) / sqrt(2.0 * M_PI * fvarDirect);
        varianceVal[it] = fvarDirect;
        
        // Computes the polar value
        kernelValPolar[it] = nik.value(t);
    }

    // Plot
    Gnuplot gp("gnuplot -persist");
    gp << "set term wxt 1\n";
    gp << "plot '-' title \"raw\" w l, '-' title \"polar\" w l\n";
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


    return 0;
}

