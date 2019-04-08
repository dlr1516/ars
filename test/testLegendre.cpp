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
#include <ars/Profiler.h>
#include <ars/ParamMap.h>
#include <boost/math/special_functions.hpp>

int main(int argc, char** argv) {
    Gnuplot gp("gnuplot -persist");
    ars::ParamMap paramMap;
    int l, m, lmax, thetaNum, plotNum, termNum;
    double theta, thetaInc;


    paramMap.read(argc, argv);
    paramMap.getParam<int>("l", l, 2);
    paramMap.getParam<int>("m", m, 1);
    paramMap.getParam<double>("theta", theta, 0.5);
    paramMap.getParam<int>("thetaNum", thetaNum, 40);
    paramMap.getParam<int>("lmax", lmax, 5);

    std::cout << "Params:\n";
    paramMap.write(std::cout);


    ars::SphericalHarmonicsLUT alegLut(lmax, thetaNum);


    plotNum = floor(2.5 * thetaNum);
    std::vector<double> thetas(plotNum + 1);
    std::vector<double> aleg1(plotNum + 1, 0.0);
    std::vector<double> aleg2(plotNum + 1, 0.0);
    double errMax = 0.0;

    thetaInc = M_PI / (plotNum);
    for (int i = 0; i <= plotNum; ++i) {
        thetas[i] = thetaInc * i;
        //xleg1[i] = ars::evaluateLegendreAssoc(l, m, xvalues[i]);
        {
            ars::ScopedTimer timer("legendreLUT");
            aleg1[i] = alegLut.evalLegendre(l, m, thetas[i]);
        }
        {
            ars::ScopedTimer timer("legendreDirect");
            aleg2[i] = sqrt(boost::math::factorial<double>(l - m) / boost::math::factorial<double>(l + m)) * boost::math::legendre_p<double>(l, m, cos(thetas[i]));
        }
        //        std::cout << "legendre at angle " << (180.0 / M_PI * thetas[i]) << " [deg], cos " << cos(thetas[i])
        //                << ": direct " << aleg2[i] << ", LUT " << aleg1[i] << ", error " << fabs(aleg1[i] - aleg2[i]) << "\n";
        if (fabs(aleg1[i] - aleg2[i]) > errMax) {
            errMax = fabs(aleg1[i] - aleg2[i]);
        }
    }
    std::cout << "\nerrMax " << errMax << "\n\nTimes (profiler):\n";
    ars::Profiler::getProfiler().printStats(std::cout);

    gp << "set term wxt " << termNum << "\n";
    gp << "plot "
            << "'-' title \"legendre LUT\" w l, "
            << "'-' title \"legendre boost\" w l"
            << "\n";
    for (int i = 0; i < thetas.size() && i < aleg1.size(); ++i) {
        gp << (180.0 / M_PI * thetas[i]) << " " << aleg1[i] << "\n";
    }
    gp << "e\n";
    for (int i = 0; i < thetas.size() && i < aleg2.size(); ++i) {
        gp << (180.0 / M_PI * thetas[i]) << " " << aleg2[i] << "\n";
    }
    gp << "e\n";
    termNum++;


    return 0;
}

