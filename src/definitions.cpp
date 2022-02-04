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

namespace ars {
    //Vec2d

    Vec2d::Vec2d() {
    }

    Vec2d::Vec2d(double v0, double v1, bool isCol = true) {
        data_[0] = v0;
        data_[1] = v1;

        isCol_ = isCol;
    }

    Vec2d::Vec2d(bool isCol = true) {
        data_[0] = 0.0;
        data_[1] = 0.0;

        isCol_ = isCol;
    }

    Vec2d::~Vec2d() {
        isCol_ = true;
    }

    void Vec2d::resetToZero() {
        data_[0] = 0.0;
        data_[1] = 0.0;
    }

    void Vec2d::multiplyByScalar(double sc) {
        data_[0] *= sc;
        data_[1] *= sc;
    }

    void Vec2d::divideByScalar(double sc) {
        data_[0] /= sc;
        data_[1] /= sc;
    }

    double Vec2d::norm() {
        return sqrt(data_[0] * data_[0] + data_[1] * data_[1]);
    }

    //VecVec2d

    VecVec2d::VecVec2d() {
    }

    VecVec2d::VecVec2d(size_t size = 0) {
        //            ptr = (cast-type*) malloc(byte-size)
        vv_ = (Vec2d*) malloc(size * sizeof (Vec2d));

        size_ = 0;
        capacity_ = size;
    }

    VecVec2d::~VecVec2d() {
        free(vv_);

        size_ = 0;
        capacity_ = 0;
    }

    void VecVec2d::pushback(Vec2d& newV) {
        if (size_ < capacity_) {
            vv_[size_] = newV;
            size_++;
        } else {
            vv_ = (Vec2d*) realloc(vv_, (capacity_ + 1) * sizeof (Vec2d));

            vv_[size_] = newV;

            size_++;
            capacity_++;
        }
    }

    void VecVec2d::pushback(const Vec2d& newV) {
        if (size_ < capacity_) {
            vv_[size_] = newV;
            size_++;
        } else {
            vv_ = (Vec2d*) realloc(vv_, (capacity_ + 1) * sizeof (Vec2d));

            vv_[size_] = newV;

            size_++;
            capacity_++;
        }
    }

    size_t VecVec2d::size() const {
        return size_;
    }

    void VecVec2d::resize(size_t sz) {
        size_ = sz;
    }

}