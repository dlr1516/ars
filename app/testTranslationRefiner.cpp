/**
 * ARS - Angular Radon Spectrum 
 * Copyright (C) 2017-2020 Dario Lodi Rizzini
 *               2021- Dario Lodi Rizzini, Ernesto Fontana
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

#include <ars/TranslationRefiner.h>

int main(int argc, char **argv) {

    std::cout << "This program executes a simple version of the standard registration algorithm Iterative Closest Point (ICP)" << std::endl;

    // define a 3 dim problem with 10000 model points
    // and 10000 template points:
    int32_t dim = 2;
    int32_t num = 100;


    ars::VectorVector2 modelPts, templatePts;

    // set model and template points
    std::cout << std::endl << "Creating model with 10000 points ..." << std::endl;
    std::cout << "Creating template by shifting model by (2,1) ..." << std::endl;
    int32_t k = 0;
    for (double x = -2; x < 2; x += 0.0004) {
        double y = 5 * x * exp(-x * x);
        ars::Vector2 ptM(x, y);
        modelPts.push_back(ptM);
        ars::Vector2 ptT(x - 2, y - 1);
        templatePts.push_back(ptT);

        k++;
        //        std::std::cout << "k " << k << std::std::endl;
    }

    // version with Eigen
    Eigen::Affine2d initialGuess;
    initialGuess.setIdentity();
    ars::TranslationRefiner translRefiner(modelPts, templatePts, initialGuess);
    Eigen::Affine2d transfOut;
    translRefiner.icpNoAssoc(transfOut);
    std::cout << "transf:" << std::endl << transfOut.matrix() << std::endl;


    return 0;
}
