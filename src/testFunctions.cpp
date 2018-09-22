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
#include <ars/functions.h>
#include <ars/ars2d.h>
#include <ars/thirdparty/gnuplot-iostream.h>

#include "ars/Profiler.h"
//#include <ars/thirdparty/Profiler.h>

struct PlotItem {
    std::vector<std::pair<double, double> > values;
    std::string title;
};

int main(int argc, char** argv) {
    std::vector<double> pnebiRecursDown;
    std::vector<double> pnebiCoeffsLUT;
    std::vector<double> pnebi0Coeffs;
    std::vector<double> arsCoeffsRecursDown;
    std::vector<double> arsCoeffsLUT;
    int n;
    double x, prec;
    double lambda, phi, cth2, sth2;

    std::cout << "Use:\n  " << argv[0] << "\n  " << argv[0] << " x\n  " << argv[0] << " x n\n" << argv[0] << " x n lut_precision\n" << std::endl;

    x = 3.0;
    n = 10;
    prec = 0.0005;
    if (argc >= 2) {
        x = atof(argv[1]);
    }
    if (argc >= 3) {
        n = atoi(argv[2]);
    }
    if (argc >= 4) {
        prec = atof(argv[3]);
    }

    std::cout << "Computing PNEBI(k," << x << ") = 2.0 * exp(-" << x << ") * besseli(k," << x << ") for k = 0, ..., " << n << std::endl;

    {
        ars::ScopedTimer timer("evaluatePNEBIVector()");
        ars::evaluatePnebiVector(n, x, pnebiRecursDown);
    }

    ars::PnebiLUT pnebiLUT(n, prec);

    std::cout << "\nPNEBI LUT summary: " << std::endl;
    pnebiLUT.printLUT(std::cout, 10, 2);
    std::cout << "-------\n\n";

    pnebiLUT.eval(x, pnebiCoeffsLUT);

    //  pnebiMatrix.resize(pnebiDirect.size(),3);
    //  pnebiMatrix << pnebiDirect, pnebiRecursDown, pnebiCoeffsLUT;
    std::cout << "PNEBI:\n";
    for (int i = 0; i <= n && i < pnebiRecursDown.size() && i < pnebiCoeffsLUT.size(); ++i) {
        std::cout << i << " \t" << pnebiRecursDown[i] << " \t" << pnebiCoeffsLUT[i] << " \n";
    }
    std::cout << std::endl;

    std::cout << "Polynom  PNEBI(0," << x << ") = " << ars::evaluatePnebi0Polynom(x) << std::endl;
    std::cout << "LUT  PNEBI(0," << x << ") = " << pnebiLUT.eval(0, x) << std::endl;

    int argnum = 1000;
    double darg = 0.05;

    std::vector<PlotItem> plotPnebiPol;
    std::vector<PlotItem> plotPnebiLUT;
    PlotItem plitem;
    int pnebiOrder = std::min(5, pnebiLUT.getOrderMax());
    for (int k = 0; k < pnebiOrder; ++k) {
        plitem.title = "pnebi_" + std::to_string(k) + "_pol";
        plotPnebiPol.push_back(plitem);
        plitem.title = "pnebi_" + std::to_string(k) + "_lut";
        plotPnebiLUT.push_back(plitem);
    }

    // Checks the value of ARS coefficients
    arsCoeffsRecursDown.resize(2 * n + 2);
    arsCoeffsLUT.resize(2 * n + 2);
    std::fill(arsCoeffsRecursDown.begin(), arsCoeffsRecursDown.end(), 0.0);
    std::fill(arsCoeffsLUT.begin(), arsCoeffsLUT.end(), 0.0);
    lambda = 2.4;
    phi = M_PI / 180.0 * (63.0);
    ars::fastCosSin(phi, cth2, sth2);
    ars::updateARSF2CoeffRecursDown(lambda, cth2, sth2, 1.0, n, arsCoeffsRecursDown);
    ars::updateARSF2CoeffRecursDownLUT(lambda, cth2, sth2, 1.0, n, pnebiLUT, arsCoeffsLUT);
    std::cout << "\nARS coefficients: lambda " << lambda << ", phi[deg] " << (180.0 / M_PI * phi) << ":" << std::endl;
    for (int i = 0; i <= 2 * n && i < arsCoeffsRecursDown.size() && i < arsCoeffsLUT.size(); ++i) {
        std::cout << i << " \t" << arsCoeffsRecursDown[i] << " \t" << arsCoeffsLUT[i] << " \n";
    }
    std::cout << std::endl;


    std::vector<double> pnebiValues(pnebiOrder);
    for (int i = 0; i < argnum; ++i) {
        double arg = darg * i;
        //    for (int k = 0; k < pnebiOrder; ++k) {
        //      assert(k < plotPnebiStd.size());
        //      plotPnebiStd[k].values.push_back( std::make_pair(arg,ans) );
        //    }
        {
            ars::ScopedTimer timer("evaluatePNEBIVector");
            ars::evaluatePnebiVector(pnebiOrder, arg, pnebiValues);
        }

        for (int k = 0; k < pnebiOrder; ++k) {
            assert(k < plotPnebiPol.size());
            plotPnebiPol[k].values.push_back(std::make_pair(arg, pnebiValues[k]));
        }
        for (int k = 0; k < pnebiOrder; ++k) {
            assert(k < plotPnebiLUT.size());
            ars::ScopedTimer timer("LUT PNEBI eval");
            double ans = pnebiLUT.eval(k, arg);
            plotPnebiLUT[k].values.push_back(std::make_pair(arg, ans));
        }
    }

    std::cout << "\nProfiler stats:\n";
    ars::Profiler::getProfiler().printStats(std::cout);
    std::cout << std::endl;

    Gnuplot gp("gnuplot -persist");
    //  std::ostream& gp = std::cout;
    double vieweps = 5e-3;
    gp << "set term wxt 0\n";
    gp << "plot ";
    for (int k = 0; k < pnebiOrder; ++k) {
        gp << "'-' title \"" << plotPnebiPol[k].title << "\" w l, "
                << "'-' title \"" << plotPnebiLUT[k].title << "\" w l";
        if (k < pnebiOrder - 1) gp << ", ";
        else gp << "\n";
    }
    for (int k = 0; k < pnebiOrder; ++k) {
        //    for (int i = 0; i < argnum; ++i) { 
        //      gp << "  " << plotPnebiStd[k].values[i].first << " " << plotPnebiStd[k].values[i].second << "\n";
        //    }
        //    gp << "e\n";
        for (int i = 0; i < argnum; ++i) {
            gp << "  " << plotPnebiPol[k].values[i].first << " " << (plotPnebiPol[k].values[i].second + vieweps) << "\n";
        }
        gp << "e\n";
        for (int i = 0; i < argnum; ++i) {
            gp << "  " << plotPnebiLUT[k].values[i].first << " " << (plotPnebiLUT[k].values[i].second) << "\n";
        }
        gp << "e\n";
    }



    //  Gnuplot gp("gnuplot -persist");
    ////  std::ostream& gp = std::cout;
    //  gp << "plot '-' title \"PNEBI0 Std.\" w l, '-' title \"PNEBI0 Polyn\" w l\n";
    //  for (int i = 0; i < argnum; ++i) { 
    //    gp << (darg * i) << " " << (funcPnebi0Standard[i]) << "\n";
    //  }
    //  gp << "e" << std::endl;
    //  for (int i = 0; i < argnum; ++i) { 
    //    gp << (darg * i) << " " << (funcPnebi0Polynom[i]) << "\n";
    //  }
    //  gp << "e" << std::endl;

    // TEST ON FAST SIN-COS FUNCTIONS
    int trigNum = 90.0;
    std::vector<double> theta(trigNum);
    std::vector<double> cosMath(trigNum);
    std::vector<double> sinMath(trigNum);
    std::vector<double> cosFast(trigNum);
    std::vector<double> sinFast(trigNum);
    double errCosMax = 0.0;
    double errSinMax = 0.0;
    for (int i = 0; i < trigNum; ++i) {
        theta[i] = 4.0 * M_PI * i / trigNum;
        cosMath[i] = cos(theta[i]);
        sinMath[i] = sin(theta[i]);
        ars::fastCosSin(theta[i], cosFast[i], sinFast[i]);
        errCosMax = std::max(errCosMax, fabs(cosMath[i] - cosFast[i]));
        errSinMax = std::max(errSinMax, fabs(sinMath[i] - sinFast[i]));
    }
    std::cout << "errCosMax " << errCosMax << ", errSinMax " << errSinMax << std::endl;
    gp << "set term wxt 1\n";
    gp << "plot '-' title \"cosMath\" w l, '-' title \"sinMath\" w l, '-' title \"cosFast\" w l, '-' title \"sinFast\" w l\n";
    for (int i = 0; i < trigNum; ++i) {
        gp << theta[i] << " " << cosMath[i] << "\n";
    }
    gp << "e\n";
    for (int i = 0; i < trigNum; ++i) {
        gp << theta[i] << " " << sinMath[i] << "\n";
    }
    gp << "e\n";
    for (int i = 0; i < trigNum; ++i) {
        gp << theta[i] << " " << cosFast[i] << "\n";
    }
    gp << "e\n";
    for (int i = 0; i < trigNum; ++i) {
        gp << theta[i] << " " << sinFast[i] << "\n";
    }
    gp << "e\n";
    return 0;
}
