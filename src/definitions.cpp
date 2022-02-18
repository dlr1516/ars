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
    //Vec2d

    //    Vec2d::Vec2d() {
    //    }
    //
    //    Vec2d::Vec2d(double v0, double v1, bool isCol = true) {
    //        data_[0] = v0;
    //        data_[1] = v1;
    //
    //        isCol_ = isCol;
    //    }
    //
    //    Vec2d::Vec2d(bool isCol = true) {
    //        data_[0] = 0.0;
    //        data_[1] = 0.0;
    //
    //        isCol_ = isCol;
    //    }
    //
    //    Vec2d::~Vec2d() {
    //        isCol_ = true;
    //    }

    //    void Vec2d::resetToZero() {
    //        data_[0] = 0.0;
    //        data_[1] = 0.0;
    //    }

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

    //    VecVec2d::VecVec2d() {
    //    }
    //
    //    VecVec2d::VecVec2d(size_t size = 0) {
    //        //            ptr = (cast-type*) malloc(byte-size)
    //        vv_ = (Vec2d*) malloc(size * sizeof (Vec2d));
    //
    //        size_ = 0;
    //        capacity_ = size;
    //    }
    //
    //    VecVec2d::~VecVec2d() {
    //        free(vv_);
    //
    //        size_ = 0;
    //        capacity_ = 0;
    //    }
    //
    //    void reserve( size_t newCapacity )
    //    {
    //        if( newCapacity < theSize )
    //            return;
    //
    //        Object *newArray = new Object[ newCapacity ];
    //        for( int k = 0; k < theSize; ++k )
    //            newArray[ k ] = std::move( objects[ k ] );
    //
    //        theCapacity = newCapacity;
    //        std::swap( objects, newArray );
    //        delete [ ] newArray;
    //    }
    //    
    //    void VecVec2d::pushback(Vec2d& newV) {
    //        if (size_ < capacity_) {
    //            vv_[size_] = newV;
    //            size_++;
    //        } else if (size_==capacity_) {
    //            vv_ = (Vec2d*) realloc(vv_, (max) * sizeof (Vec2d));
    //
    //            vv_[size_] = newV;
    //
    //            size_++;
    //            capacity_ = max;
    //        }
    //        else {
    //            ARS_ERROR("VecVec2d pushback is failing");
    //        }
    //    }
    //
    //    void VecVec2d::pushback(const Vec2d& newV) {
    //        if (size_ < capacity_) {
    //            vv_[size_] = newV;
    //            size_++;
    //        } else if (size_==capacity_) {
    //            size_t max = 1 > (capacity_*2) ? 1 : (capacity_*2);  
    //            vv_ = (Vec2d*) realloc(vv_, (max) * sizeof (Vec2d));
    //
    //            vv_[size_] = newV;
    //
    //            size_++;
    //            capacity_ = max;
    //        }
    //        else {
    //            ARS_ERROR("VecVec2d pushback is failing");
    //        }
    //    }
    //
    //    size_t VecVec2d::size() const {
    //        return size_;
    //    }
    //
    //    void VecVec2d::resize(size_t sz) {
    //        size_ = sz;
    //    }
    //

    //Mat2d

    //    Mat2d::Mat2d() {
    //        data_[0] = 0.0;
    //        data_[1] = 0.0;
    //        data_[2] = 0.0;
    //        data_[3] = 0.0;
    //    }

    //    Mat2d::~Mat2d() {
    //    }

    //    void Mat2d::resetToZero() {
    //        data_[0 * Two + 0] = 0.0; // = data[0]
    //        data_[0 * Two + 1] = 0.0; // = data[1]
    //        data_[1 * Two + 0] = 0.0; // = data[2]
    //        data_[1 * Two + 1] = 0.0; // = data[3]
    //    }

    void Mat2d::setToIdentity() {
        data_[0 * Two + 0] = 1.0; // = data[0]
        data_[0 * Two + 1] = 0.0; // = data[1]
        data_[1 * Two + 0] = 0.0; // = data[2]
        data_[1 * Two + 1] = 1.0; // = data[3]
    }

    void Mat2d::multiplyByScalar(double sc) {
        data_[0] *= sc;
        data_[1] *= sc;
        data_[2] *= sc;
        data_[3] *= sc;
    }

    void Mat2d::divideByScalar(double sc) {
        data_[0] /= sc;
        data_[1] /= sc;
        data_[2] /= sc;
        data_[3] /= sc;
    }

    //    void Mat2d::fillRowMajor(double a, double b, double c, double d) {
    //        data_[0 * Two + 0] = a;
    //        data_[0 * Two + 1] = b;
    //        data_[1 * Two + 0] = c;
    //        data_[1 * Two + 1] = d;
    //    }

    //    void Mat2d::make2dRotMat(double theta) {
    //        data_[0 * Two + 0] = cos(theta);
    //        data_[0 * Two + 1] = -sin(theta); //avoiding useless function calling
    //        data_[1 * Two + 0] = -data_[0 * Two + 1];
    //        data_[1 * Two + 1] = data_[0 * Two + 0];
    //    }

    //    void Mat2d::transpose() {
    //        double tmp = data_[0 * Two + 1];
    //        data_[0 * Two + 1] = data_[1 * Two + 0];
    //        data_[1 * Two + 0] = tmp;
    //    }

    Mat2d Mat2d::transposeReturningValue() {
        Mat2d transposed;
        transposed.data_[0 * Two + 0] = data_[0 * Two + 0];
        transposed.data_[0 * Two + 1] = data_[0 * Two + 1];
        transposed.data_[1 * Two + 0] = data_[1 * Two + 0];
        transposed.data_[1 * Two + 1] = data_[1 * Two + 1];

        return transposed;
    }

    double Mat2d::determinant() const {
        return data_[0 * Two + 0] * data_[1 * Two + 1] - data_[0 * Two + 1] * data_[1 * Two + 0];
    }

    double Mat2d::trace() const {
        return data_[0 * Two + 0] + data_[1 * Two + 1];
    }

    void Mat2d::invert() {
        //            double det = data_[0 * Two + 0] * data_[1 * Two + 1] - data_[0 * Two + 1] * data_[1 * Two + 0]; //maybe use directly determinant() function??
        double detInv = 1.0 / determinant();

        double aOrig = data_[0 * Two + 0];
        double bOrig = data_[0 * Two + 1];
        double cOrig = data_[1 * Two + 0];
        double dOrig = data_[1 * Two + 1];

        data_[0 * Two + 0] = dOrig * detInv;
        data_[0 * Two + 1] = -bOrig * detInv;
        data_[1 * Two + 0] = -cOrig * detInv;
        data_[1 * Two + 1] = aOrig * detInv;

    }

    Mat2d Mat2d::inverse() const {
        Mat2d m;

        double detInv = 1.0 / determinant();

        double aOrig = data_[0 * Two + 0];
        double bOrig = data_[0 * Two + 1];
        double cOrig = data_[1 * Two + 0];
        double dOrig = data_[1 * Two + 1];

        m.data_[0 * Two + 0] = dOrig * detInv;
        m.data_[0 * Two + 1] = -bOrig * detInv;
        m.data_[1 * Two + 0] = -cOrig * detInv;
        m.data_[1 * Two + 1] = aOrig * detInv;

        return m;
    }

    //    void Mat2d::setDiagonal(double a11, double a22) {
    //        data_[0 * Two + 0] = a11;
    //        data_[1 * Two + 1] = a22;
    //    }

}