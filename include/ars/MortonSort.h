#ifndef MORTONSORT_H
#define MORTONSORT_H

#include <cmath>
#include <cstdint>
#include <Eigen/Dense>

namespace ars {
    
    // ---------------------------------------------------------------
    // NUMBER OF LEADING ZEROS
    // ---------------------------------------------------------------
    
    int8_t nlz8(uint8_t x);

    int16_t nlz16(uint16_t x);

    int32_t nlz32(uint32_t x);

    int64_t nlz64(uint64_t x);
    
    // ---------------------------------------------------------------
    // MORTON INTEGER TRAITS (FOR TREATING THE SIGN OF INTEGER)
    // ---------------------------------------------------------------

    template <typename I> struct MortonIntegerTraits;

    template <> struct MortonIntegerTraits<int8_t> {
        using IntegerType = int8_t;
        using UnsignedType = uint8_t;

        static const int BIT_NUM = 8;

        static UnsignedType removeSign(const IntegerType& i) {
            return (i ^ 0x01);
        }
        
        static IntegerType nlz(const IntegerType& i) {
            return nlz8(i);
        }
    };

    template <> struct MortonIntegerTraits<uint8_t> {
        using IntegerType = uint8_t;
        using UnsignedType = uint8_t;

        static const int BIT_NUM = 8;

        static UnsignedType removeSign(const IntegerType& i) {
            return i;
        }
        
        static IntegerType nlz(const IntegerType& i) {
            return nlz8(i);
        }
    };

    template <> struct MortonIntegerTraits<int16_t> {
        using IntegerType = int16_t;
        using UnsignedType = uint16_t;

        static const int BIT_NUM = 16;

        static UnsignedType removeSign(const IntegerType& i) {
            return (i ^ 0x0001);
        }
        
        static IntegerType nlz(const IntegerType& i) {
            return nlz16(i);
        }
    };

    template <> struct MortonIntegerTraits<uint16_t> {
        using IntegerType = uint16_t;
        using UnsignedType = uint16_t;

        static const int BIT_NUM = 16;

        static UnsignedType removeSign(const IntegerType& i) {
            return i;
        }
        
        static IntegerType nlz(const IntegerType& i) {
            return nlz16(i);
        }
    };

    template <> struct MortonIntegerTraits<int32_t> {
        using IntegerType = int32_t;
        using UnsignedType = uint32_t;

        static const int BIT_NUM = 32;

        static UnsignedType removeSign(const IntegerType& i) {
            return (i ^ 0x00000001);
        }
        
        static IntegerType nlz(const IntegerType& i) {
            return nlz32(i);
        }
    };

    template <> struct MortonIntegerTraits<uint32_t> {
        using IntegerType = uint32_t;
        using UnsignedType = uint32_t;

        static const int BIT_NUM = 32;

        static UnsignedType removeSign(const IntegerType& i) {
            return i;
        }
        
        static IntegerType nlz(const IntegerType& i) {
            return nlz32(i);
        }
    };

    template <> struct MortonIntegerTraits<int64_t> {
        using IntegerType = int64_t;
        using UnsignedType = uint64_t;

        static const int BIT_NUM = 64;

        static UnsignedType removeSign(const IntegerType& i) {
            return (i ^ 0x0000000000000001);
        }
        
        static IntegerType nlz(const IntegerType& i) {
            return nlz64(i);
        }
    };

    template <> struct MortonIntegerTraits<uint64_t> {
        using IntegerType = uint64_t;
        using UnsignedType = uint64_t;

        static const int BIT_NUM = 64;

        static UnsignedType removeSign(const IntegerType& i) {
            return i;
        }
        
        static IntegerType nlz(const IntegerType& i) {
            return nlz64(i);
        }
    };

    // ---------------------------------------------------------------
    // MORTON COMPARATOR
    // ---------------------------------------------------------------

    /**
     * MSB: Most Significant Bit comparison. 
     * Given two integers i1 and i2, it computes if floor(log2(i1)) < floor(log2(i2))
     * using the well known Timothy Chan's trick illustrated in 
     *   
     * T.M. Chan, "A minimalist's implementation of an approximate nearest 
     * neighbor algorithm in fixed dimensions", 2006
     * http://tmc.web.engr.illinois.edu/pub_ann.html
     * 
     * E.g. i1 = 8 = 1000b, i2 = 11 = 1011b, i1 ^ i2 = 0011
     *    i1 < i2 -> TRUE, i1 = 1000b < 
     * 
     * @param i1 the first integer
     * @param i2 the second integer
     * @return 
     */
    template <typename Integer>
    inline bool msb(const Integer& i1, const Integer& i2) {
        return (i1 < i2 && i1 < (i1 ^ i2));
    }

    /**
     * Compares two (Eigen) vectors of integer values with integer type I and
     * dimension Dim and sort them according to Morton order.
     * It is based on the so called Chan's trick:
     * 
     * T.M. Chan, "A minimalist's implementation of an approximate nearest 
     * neighbor algorithm in fixed dimensions", 2006
     * http://tmc.web.engr.illinois.edu/pub_ann.html
     * 
     * @param v1 the first integer vector
     * @param v2 the second integer vector
     * @return true if v1 is before v2 in Morton Order. 
     */
    template <typename I, int Dim>
    bool mortonCmpInt(const Eigen::Matrix<I, Dim, 1>& v1, const Eigen::Matrix<I, Dim, 1>& v2) {
        using MIT = MortonIntegerTraits<I>;
        using UnsignedType = typename MIT::UnsignedType;
        UnsignedType lastDim, lastXor, currXor;
        lastDim = 0;
        lastXor = MIT::removeSign(v1(0)) ^ MIT::removeSign(v2(0));
        for (int d = 1; d < Dim; ++d) {
            currXor = MIT::removeSign(v1(d)) ^ MIT::removeSign(v2(d));
            if (msb<UnsignedType>(lastXor, currXor)) {
                lastDim = d;
                lastXor = currXor;
            }
        }
        return (v1(lastDim) < v2(lastDim));
    }
    
    template <typename I, int Dim>
    int mortonLevel(const Eigen::Matrix<I, Dim, 1>& v1, const Eigen::Matrix<I, Dim, 1>& v2) {
        using MIT = MortonIntegerTraits<I>;
        using UnsignedType = typename MIT::UnsignedType;
        typename MIT::IntegerType levelMax, level;
        levelMax = 0;
        for (int d = 0; d < Dim; ++d) {
            level = MIT::BIT_NUM - MIT::nlz(MIT::removeSign(v1(d)) ^ MIT::removeSign(v2(d)));
            if (level > levelMax) {
                levelMax = level;
            }
        }
        return levelMax;
    }

}

#endif /* MORTONSORT_H */

