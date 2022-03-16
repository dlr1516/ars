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

#include <ars/definitions.h>

namespace cuars {

    Affine2d::Affine2d() {
        rot_ = 0.0;
        translX_ = 0.0;
        translY_ = 0.0;

        initdata(rot_, translX_, translY_);
    }

    Affine2d::Affine2d(double rot, double tx, double ty) {
        rot_ = rot;
        translX_ = tx;
        translY_ = ty;

        initdata(rot, tx, ty);
    }

    Affine2d::~Affine2d() {
    }

    void Affine2d::initdata(double r, double tx, double ty) {
        data_[0 * Three + 0] = cos(r);
        data_[0 * Three + 1] = -sin(r);
        data_[1 * Three + 0] = -data_[0 * Three + 1];
        data_[1 * Three + 1] = data_[0 * Three + 0];

        data_[0 * Three + 2] = tx;
        data_[1 * Three + 2] = ty;

        data_[2 * Three + 2] = 1.0;

        data_[2 * Three + 0] = 0.0;
        data_[2 * Three + 1] = 0.0;
    }

    bool Affine2d::isLastRowOK() const {
        double a20 = data_[2 * Three + 0];
        double a21 = data_[2 * Three + 1];
        double a22 = data_[2 * Three + 2];

        if (a20 == 0 && a21 == 0 && a22 == 1)
            return true;

        printf("BAD LAST ROW\n");
        return false;
    }

    bool Affine2d::isScale1() {
        double a22 = data_[2 * Three + 2];

        if (a22 == 1)
            return true;

        printf("BAD SCALE\n");
        return false;
    }

    double Affine2d::at(int r, int c) const {
        if (r >= 0 && r < Three && c >= 0 && c < Three)
            return data_[r * Three + c];
        else {
            printf("ERROR accessing matrix with .at() method!\n");
            return 1000000;
        }

    }

}