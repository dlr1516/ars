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
#include <Eigen/Dense>
#include <iostream>

#include <ars/FourierLowerUpperLut.h>
#include <ars/ParamMap.h>
#include <ars/definitions.h>
#include <ars/functions.h>

#include <ars/thirdparty/gnuplot-iostream.h>
#include <chrono>

#define PRINT_DIM(X)                                                \
    std::cout << #X << " rows " << X.rows() << " cols " << X.cols() \
              << std::endl;
#define RAD2DEG(X) (180.0 / M_PI * (X))

using BoundInterval = ars::FourierLowerUpperLut::Interval;

// std::vector<double> coeffs = {0.0956209,    0,
//                               -0.00756418,  0.0254858,
//                               -0.0159511,   -0.0222764,
//                               0.0105998,    -0.00548222,
//                               -0.00154664,  0.0105675,
//                               0.00202689,   0.00236173,
//                               0.00606581,   0.00169872,
//                               0.000362685,  -0.000263201,
//                               0.00168924,   -0.00354879,
//                               0.000656388,  -0.00223724,
//                               -0.00162124,  -9.29845e-06,
//                               -0.00127703,  -0.0032378,
//                               0.00142952,   0.00148903,
//                               -0.000444384, 0.000975841,
//                               0.000341318,  -0.00198222,
//                               1.51505e-05,  0.000739245,
//                               -0.000620007, 0.000517741,
//                               -0.000155,    -0.000509777,
//                               0.000594383,  -0.000966466,
//                               -0.000169302, -0.000123654,
//                               -0.000989125, 0.000505298};

std::vector<double> coeffs = {0.0956209,  0,          -0.00756418, 0.0254858,
                              -0.0159511, -0.0222764, 0.0105998,   -0.00548222};

void plotBranchBoundBox(std::ostream& out,
                        const std::vector<BoundInterval>& bbbs);

int main(int argc, char** argv) {
    ars::FourierLowerUpperLut lut;
    ars::ParamMap params;
    int levelNum;
    int thetaNum;
    double xL, xU;

    params.read(argc, argv);
    params.getParam<int>("level_num", levelNum, 4);
    params.getParam<int>("theta_num", thetaNum, 180);
    params.getParam<double>("theta_l", xL, 30.0);
    params.getParam<double>("theta_u", xU, 50.0);
    xL *= (M_PI / 180.0);
    xU *= (M_PI / 180.0);

    std::cout << "Parameters: " << std::endl;
    params.write(std::cout);

    lut.init(coeffs, levelNum);

    // double yL, yU;
    // lut.findLU(xL, xU, yL, yU);
    // std::cout << "Bounds for [" << xL << ", " << xU << "]: " << yL << ", " <<
    // yU
    //           << std::endl;

    Gnuplot gp("gnuplot -persist");
    double vieweps = 5e-3;

    lut.exportPlot(gp);

    //  std::ostream& gp = std::cout;
    gp << "set term wxt 100\n";
    gp << "plot '-' title \"fourier\" w l, '-' title "
          "\"bb\" w l\n";
    for (int i = 0; i < thetaNum; ++i) {
        double theta = (M_PI / thetaNum) * i;
        double fourier = ars::evaluateFourier(coeffs, 2.0 * theta);
        gp << RAD2DEG(theta) << " " << fourier << "\n";
    }
    gp << "e" << std::endl;
    plotBranchBoundBox(gp, lut.intervals());
    gp << "e" << std::endl;

    return 0;
}

void plotBranchBoundBox(std::ostream& out,
                        const std::vector<BoundInterval>& bbbs) {
    for (auto& bbb : bbbs) {
        out << RAD2DEG(bbb.thetaMin) << " " << bbb.yLower << "\n"
            << RAD2DEG(bbb.thetaMax) << " " << bbb.yLower << "\n"
            << RAD2DEG(bbb.thetaMax) << " " << bbb.yUpper << "\n"
            << RAD2DEG(bbb.thetaMin) << " " << bbb.yUpper << "\n"
            << RAD2DEG(bbb.thetaMin) << " " << bbb.yLower << "\n\n";
    }
}
