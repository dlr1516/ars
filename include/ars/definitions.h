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
#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <iostream>
#include <cmath>

#define ARS_PRINT(MSG) std::cout << __FILE__ << "," << __LINE__ << ": " << MSG << std::endl;

#define ARS_ERROR(MSG) std::cerr << __FILE__ << "," << __LINE__ << ": " << MSG << std::endl;

#define ARS_VARIABLE(X1) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) << std::endl;

#define ARS_VARIABLE2(X1,X2) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) \
   << ", " << (#X2) << " " << (X2) << std::endl;

#define ARS_VARIABLE3(X1,X2,X3) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) \
   << ", " << (#X2) << " " << (X2) << ", " << (#X3) << " " << (X3) << std::endl;

#define ARS_VARIABLE4(X1,X2,X3,X4) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) \
   << ", " << (#X2) << " " << (X2) << ", " << (#X3) << " " << (X3) << ", " << (#X4) << " " << (X4) << std::endl;

#define ARS_VARIABLE5(X1,X2,X3,X4,X5) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) \
   << ", " << (#X2) << " " << (X2) << ", " << (#X3) << " " << (X3) << ", " << (#X4) << " " << (X4) \
   << ", " << (#X5) << " " << (X5) << std::endl;

#define ARS_VARIABLE6(X1,X2,X3,X4,X5,X6) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) \
   << ", " << (#X2) << " " << (X2) << ", " << (#X3) << " " << (X3) << ", " << (#X4) << " " << (X4) \
   << ", " << (#X5) << " " << (X5) << ", " << (#X6) << " " << (X6) << std::endl;

#define ARS_VARIABLE7(X1,X2,X3,X4,X5,X6,X7) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) \
   << ", " << (#X2) << " " << (X2) << ", " << (#X3) << " " << (X3) << ", " << (#X4) << " " << (X4) \
   << ", " << (#X5) << " " << (X5) << ", " << (#X6) << " " << (X6) << ", " << (#X7) << " " << (X7) << std::endl;

#define ARS_ASSERT(COND) if (!(COND)) { std::cerr << __FILE__ << "," << __LINE__ << ": assertion failed on " << #COND << std::endl; exit(-1); }


namespace ars {

    static const size_t Two = 2; //useful for expanding (i,j) indexing into  i*Two+j

    //    using Vector2 = Eigen::Vector2d;

    class Vec2d {
    public:
        double data_[2];

        void resetToZero() {
            data_[0] = 0.0;
            data_[1] = 0.0;
        }

        double norm() {
            return sqrt(data_[0] * data_[0] + data_[1] * data_[1]);
        }
    };

    //    using VectorVector2 = std::vector<Vector2, Eigen::aligned_allocator<Vector2> >    

    class VecVec2d {
    private:
        size_t size_; //number of elements
    public:
        Vec2d *vv_;

        size_t size() const {
            return size_;
        }
    };



    //    using Matrix2 = Eigen::Matrix2d;

    class Mat2d {
    public:
        double data_[4];
        
        void resetToZero() {
            data_[0 * Two + 0] = 0.0; // = data[0]
            data_[0 * Two + 1] = 0.0; // = data[1]
            data_[1 * Two + 0] = 0.0; // = data[2]
            data_[1 * Two + 1] = 0.0; // = data[3]
        }

        void fillRowMajor(double a, double b, double c, double d) {
            data_[0 * Two + 0] = a;
            data_[0 * Two + 1] = b;
            data_[1 * Two + 0] = c;
            data_[1 * Two + 1] = d;
        }

        void make2dRotMat(double theta) {
            data_[0 * Two + 0] = cos(theta);
            data_[0 * Two + 1] = -sin(theta); //avoiding useless function calling
            data_[1 * Two + 0] = -data_[0 * Two + 1];
            data_[1 * Two + 1] = data_[0 * Two + 0];
        }

        void transpose() {
            double tmp = data_[0 * Two + 1];
            data_[0 * Two + 1] = data_[1 * Two + 0];
            data_[1 * Two + 0] = tmp;
        }

        void setDiagonal(double a11, double a22) {
            data_[0 * Two + 0] = a11;
            data_[1 * Two + 1] = a22;
        }       
    };

    //    using VectorMatrix2 = std::vector<Matrix2, Eigen::aligned_allocator<Matrix2> >;

    class VecMat2d {
    private:
        size_t size_; //number of elements
    public:
        Mat2d *mm_;

        size_t size() const {
            return size_;
        }
    };




} // end of namespace

#endif /* DEFINITIONS_H */

