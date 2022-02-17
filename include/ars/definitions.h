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
#include <limits>
#include <thrust/host_vector.h>
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


namespace cuars {

    static const size_t Two = 2; //useful for expanding (i,j) indexing into  i*Two+j
    static const size_t Three = 3; //useful for expanding (i,j) indexing into  i*Two+j
    static const size_t Nine = 9; //useful for expanding (i,j) indexing into  i*Two+j



    //    using Vector2 = Eigen::Vector2d;

    using Vec2d = double2;

    //    class Vec2d {
    //    public:
    //        double data_[2];
    //        bool isCol_; //default -> true
    //
    //        //        Vec2d() {
    //        //            data_[0] = 0.0;
    //        //            data_[1] = 0.0;
    //        //
    //        //            isCol_ = true;
    //        //        }
    //
    //        Vec2d();
    //
    //        Vec2d(double v0, double v1, bool isCol);
    //
    //        Vec2d(bool isCol);
    //
    //        virtual ~Vec2d();
    //
    //        void resetToZero();
    //
    //        void multiplyByScalar(double sc);
    //
    //        void divideByScalar(double sc);
    //
    //        double norm();
    //    };

    using VecVec2d = thrust::host_vector<Vec2d>;

    //    class VecVec2d {
    //    private:
    //        size_t size_; //number of elements
    //        size_t capacity_;
    //    public:
    //        Vec2d *vv_;
    //
    //        VecVec2d();
    //        
    //        VecVec2d(size_t size);
    //
    //        virtual ~VecVec2d();
    //
    //        void pushback(Vec2d& newV);
    //
    //        void pushback(const Vec2d& newV);
    //
    //        size_t size() const;
    //
    //        void resize(size_t sz);
    //    };



    //    using Matrix2 = Eigen::Matrix2d;

    //    class Mat2d {
    //    public:
    //        double data_[4];
    //
    //        Mat2d();
    //
    //        virtual ~Mat2d();
    //
    //        void resetToZero();
    //
    //        void setToIdentity();
    //
    //        void multiplyByScalar(double sc);
    //
    //        void divideByScalar(double sc);
    //
    //        void fillRowMajor(double a, double b, double c, double d);
    //
    //        void make2dRotMat(double theta);
    //
    //        void transpose();
    //
    //        Mat2d transposeReturningValue();
    //
    //        double determinant() const;
    //
    //        double trace() const;
    //
    //        void invert();
    //
    //        Mat2d inverse() const;
    //
    //        void setDiagonal(double a11, double a22);
    //
    //    };

    using Mat2d = double4; //w x \n y z

    using VecMat2d = thrust::host_vector<Mat2d>;

    //    class VecMat2d {
    //    private:
    //        size_t size_; //number of elements
    //        size_t capacity_;
    //
    //    public:
    //        Mat2d *mm_;
    //
    //        VecMat2d(size_t size = 0) {
    //            //            ptr = (cast-type*) malloc(byte-size)
    //            mm_ = (Mat2d*) malloc(size * sizeof (Mat2d));
    //
    //            size_ = 0;
    //            capacity_ = size;
    //        }
    //
    //        virtual ~VecMat2d() {
    //            free(mm_);
    //
    //            size_ = 0;
    //            capacity_ = 0;
    //        }
    //
    //        void pushback(Mat2d& newM) {
    //            if (size_ < capacity_) {
    //                mm_[size_] = newM;
    //                size_++;
    //            } else {
    //                mm_ = (Mat2d*) realloc(mm_, (capacity_ + 1) * sizeof (Vec2d));
    //
    //                mm_[size_] = newM;
    //
    //                size_++;
    //                capacity_++;
    //            }
    //        }
    //
    //        void pushback(const Mat2d& newM) {
    //            if (size_ < capacity_) {
    //                mm_[size_] = newM;
    //                size_++;
    //            } else {
    //                mm_ = (Mat2d*) realloc(mm_, (capacity_ + 1) * sizeof (Vec2d));
    //
    //                mm_[size_] = newM;
    //
    //                size_++;
    //                capacity_++;
    //            }
    //        }
    //
    //        size_t size() const {
    //            return size_;
    //        }
    //
    //        void resize(size_t sz) {
    //            size_ = sz;
    //        }
    //    };

    struct Affine2d {
        double data_[Nine];

        double rot_;
        double translX_;
        double translY_;

        Affine2d() {
        }

        Affine2d(double rot, double tx, double ty) {
            rot_ = rot;
            translX_ = tx;
            translY_ = ty;

            initdata(rot, tx, ty);

        }

        void initdata(double r, double tx, double ty) {
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

        bool isLastRowOK() {
            double a20 = data_[2 * Three + 0];
            double a21 = data_[2 * Three + 1];
            double a22 = data_[2 * Three + 2];

            if (a20 == 0 && a21 == 0 && a22 == 1)
                return true;

            printf("BAD LAST ROW\n");
            return false;
        }

        bool isScale1() {
            double a22 = data_[2 * Three + 2];

            if (a22 == 1)
                return true;

            printf("BAD SCALE\n");
            return false;
        }

        double at(int r, int c) {
            if (r >= 0 && r < Three && c >= 0 && c < Three)
                return data_[r * Three + c];
            else {
                printf("ERROR accessing matrix with .at() method!\n");
                return 1000000;
            }

        }

    };

} // end of namespace

#endif /* DEFINITIONS_H */

