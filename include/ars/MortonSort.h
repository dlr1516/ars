#ifndef MORTONSORT_H
#define MORTONSORT_H

#include <cmath>
#include <stdint>

namespace ars {

    void getExpMantissaF(const float& f, uint32_t& e, uint32_t& m, bool& s) {

        union {
            float f;
            uint32_t i;
        } floatBits;

        floatBits.f = f;
        m = 0x007FFFFF & floatBits.i;
        e = 0x7F800000 & floatBits.i;
        s = 0x80000000 & floatBits.i;
    }

    void getExpMantissaD(const double& f, uint64_t& e, uint64_t& m, bool& s) {

        union {
            double f;
            uint64_t i;
        } floatBits;

        floatBits.f = f;
        m = 0x000FFFFFFFFFFFFF & floatBits.i;
        e = 0x7FF0000000000000 & floatBits.i;
        s = 0x8000000000000000 & floatBits.i;
    }

    // --------------------------------------------------------------
    // XOR MSB 
    // --------------------------------------------------------------

    template <typename Float> struct XorMsb;

    template <>
    struct XorMsb<float> {
        using Float = float;
        using Integer = int;

        union FloatInteger {
            Float f;
            Integer i;
        };

        static Integer eval(const Float& f1, const Float& f2) {
            Integer e1, e2, ret;
            bool s1, s2;
            FloatInteger u1, u2, zero;

            //getExpMantissaF(f1, m1, e1, s1);
            //getExpMantissaF(f2, m2, e2, s2);
            u1.f = frexpf(f1, &e1);
            u2.f = frexpf(f2, &e2);
            zero.f = 0.0f;

            if (e1 == e2) {
                if (m1 == m2) {
                    ret = 0;
                } else {
                    u1.i = ((u1.i ^ u2.i) | zero.i);
                    u1.f = frexpf(u1.f - 0.5f, &e1);
                    ret = e1 + e2;
                }
            } else {
                if (e1 < e2)
                    ret = e2;
                else
                    ret = e1;
            }
            return ret;
        }
    };

    /**
     * Given two integers i1 and i2, it computes if floor(log2(i1)) < floor(log2(i2))
     * using the well known Timothy Chan's trick illustrated in 
     *   
     * T.M. Chan, "A minimalist's implementation of an approximate nearest 
     * neighbor algorithm in fixed dimensions", 2006
     * http://tmc.web.engr.illinois.edu/pub_ann.html
     * 
     * @param i1 the first integer
     * @param i2 the second integer
     * @return 
     */
    template <typename Integer>
    bool lessMbs(const Integer& i1, const Integer& i2) {
        return (i1 < i2 && i1 < (i1 ^ i2));
    }

    bool xorMsbD(const double& f1, const double& f2) {
        int e1, e2;
        double m1, m2;
        m1 = frexp(f1, e1);
        m2 = frexp(f2, e2);
        if (e1 == e2) {

        }
    }

    template <typename Integer, typename Dim>
    bool mortonCmpInt(const Eigen::Matrix<Integer, Dim, 1>& v1, const Eigen::Matrix<Integer, Dim, 1>& v2) {
        Integer dimLast = 0;
        Integer diffLast = 0;
        Integer diffShifted;
        for (Integer d = 0; d < Dim; ++d) {
            diffShifted = v1(d) ^ v2(d);
            if (lessMbs(diffLast, diffShifted)) {
                dimLast = d;
                diffLast = diffShifted;
            }
        }
        return v1(dimLast) < v2(dimLast);
    }

    template <typename Float, typename Dim>
    bool mortonCmpFloat(const Eigen::Matrix<Float, Dim, 1>& v1, const Eigen::Matrix<Float, Dim, 1>& v2) {
        Float x, y;
        int dim;

        dim = 0;
        x = 0;
        for (int d = 0; d < Dim; ++d) {
            y = xorMsb(v1(d), v2(d));
            if (x < y) {
                x = y;
                dim = d;
            }
        }
        return v1(dim) < v2(dim);
    }

    template <typename Scalar, typename Dim>
    struct MortonSort {
        using Point = Eigen::Matrix<Scalar, Dim, 1>;

        static bool less(const Point& p1, const Point& p2);
    };

    template <typename Dim>
    struct MortonSort<int, Dim> {
        using Point = Eigen::Matrix<int, Dim, 1>;

        static bool less(const Point& p1, const Point& p2) {

        }
    };

}

#endif /* MORTONSORT_H */

