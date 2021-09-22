#ifndef MORTONSORT_H
#define MORTONSORT_H

#include <bitset>
#include <cmath>
#include <cstdint>

#include "definitions.h"

namespace ars {

    // ---------------------------------------------------------------
    // NUMBER OF LEADING ZEROS
    // ---------------------------------------------------------------

    int8_t nlz8(uint8_t x);

    int16_t nlz16(uint16_t x);

    int32_t nlz32(uint32_t x);

    int64_t nlz64(uint64_t x);
    
    // ---------------------------------------------------------------
    // FLOOR/CEIL LOW POWER 2
    // ---------------------------------------------------------------
    
    /**
     * Returns the larger power of 2 less than the given argument.
     * Example: flp2(5) = 2^2.
     * @param x input argument
     */
    uint8_t flp2u8(uint8_t x);

    /**
     * Returns the larger power of 2 less than the given argument.
     * Example: flp2(5) = 2^2.
     * @param x input argument
     */
    uint16_t flp2u16(uint16_t x);

    /**
     * Returns the larger power of 2 less than the given argument.
     * Example: flp2(5) = 2^2.
     * @param x input argument
     */
    uint32_t flp2u32(uint32_t x);

    /**
     * Returns the larger power of 2 less than the given argument.
     * Example: flp2(5) = 2^2.
     * @param x input argument
     */
    uint64_t flp2u64(uint64_t x);

    /**
     * Returns the smaller power of 2 greater than the given argument.
     * Example: clp2(5) = 2^3.
     * @param x input argument
     */
    uint8_t clp2u8(uint8_t x);

    /**
     * Returns the smaller power of 2 greater than the given argument.
     * Example: clp2(5) = 2^3.
     * @param x input argument
     */
    uint16_t clp2u16(uint16_t x);

    /**
     * Returns the smaller power of 2 greater than the given argument.
     * Example: clp2(5) = 2^3.
     * @param x input argument
     */
    uint32_t clp2u32(uint32_t x);

    /**
     * Returns the smaller power of 2 greater than the given argument.
     * Example: clp2(5) = 2^3.
     * @param x input argument
     */
    uint64_t clp2u64(uint64_t x);

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

        static IntegerType log2Mod(const IntegerType& i) {
            if (i > 0)
                return (BIT_NUM - 1 - nlz8(i));
            else if (i < 0)
                return (BIT_NUM - 1 - nlz8(-i));
            else
                return 0; // arbitrary output for invalid argument 0
        }

        static IntegerType floorPow2(const IntegerType& i) {
            if (i >= 0)
                return flp2u8((UnsignedType) i);
            else
                return -clp2u8((UnsignedType) (-i));
        }

        static IntegerType ceilPow2(const IntegerType& i) {
            if (i >= 0)
                return clp2u8((UnsignedType) i);
            else
                return -flp2u8((UnsignedType) (-i));
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

        static IntegerType log2Mod(const IntegerType& i) {
            if (i != 0)
                return (BIT_NUM - nlz8(i) - 1);
            else
                return 0; // arbitrary output for invalid argument 0
        }

        static IntegerType floorPow2(const IntegerType& i) {
            return flp2u8(i);
        }

        static IntegerType ceilPow2(const IntegerType& i) {
            return clp2u8(i);
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

        static IntegerType log2Mod(const IntegerType& i) {
            if (i > 0)
                return (BIT_NUM - 1 - nlz16(i));
            else if (i < 0)
                return (BIT_NUM - 1 - nlz16(-i));
            else
                return 0; // arbitrary output for invalid argument 0
        }

        static IntegerType floorPow2(const IntegerType& i) {
            if (i >= 0)
                return flp2u16((UnsignedType) i);
            else
                return -clp2u16((UnsignedType) (-i));
        }

        static IntegerType ceilPow2(const IntegerType& i) {
            if (i >= 0)
                return clp2u16((UnsignedType) i);
            else
                return -flp2u16((UnsignedType) (-i));
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

        static IntegerType log2Mod(const IntegerType& i) {
            if (i != 0)
                return (BIT_NUM - nlz16(i) - 1);
            else
                return 0; // arbitrary output for invalid argument 0
        }

        static IntegerType floorPow2(const IntegerType& i) {
            return flp2u16(i);
        }

        static IntegerType ceilPow2(const IntegerType& i) {
            return clp2u16(i);
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

        static IntegerType log2Mod(const IntegerType& i) {
            if (i > 0)
                return (BIT_NUM - 1 - nlz32(i));
            else if (i < 0)
                return (BIT_NUM - 1 - nlz32(-i));
            else
                return 0; // arbitrary output for invalid argument 0
        }

        static IntegerType floorPow2(const IntegerType& i) {
            if (i >= 0)
                return flp2u32((UnsignedType) i);
            else
                return -clp2u32((UnsignedType) (-i));
        }

        static IntegerType ceilPow2(const IntegerType& i) {
            if (i >= 0)
                return clp2u32((UnsignedType) i);
            else
                return -flp2u32((UnsignedType) (-i));
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

        static IntegerType log2Mod(const IntegerType& i) {
            if (i != 0)
                return (BIT_NUM - nlz32(i) - 1);
            else
                return 0; // arbitrary output for invalid argument 0
        }

        static IntegerType floorPow2(const IntegerType& i) {
            return flp2u32(i);
        }

        static IntegerType ceilPow2(const IntegerType& i) {
            return clp2u32(i);
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

        static IntegerType log2Mod(const IntegerType& i) {
            if (i > 0)
                return (BIT_NUM - 1 - nlz8(i));
            else if (i < 0)
                return (BIT_NUM - 1 - nlz8(-i));
            else
                return 0; // arbitrary output for invalid argument 0
        }

        static IntegerType floorPow2(const IntegerType& i) {
            if (i >= 0)
                return flp2u64((UnsignedType) i);
            else
                return -clp2u64((UnsignedType) (-i));
        }

        static IntegerType ceilPow2(const IntegerType& i) {
            if (i >= 0)
                return clp2u64((UnsignedType) i);
            else
                return -flp2u64((UnsignedType) (-i));
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

        static IntegerType log2Mod(const IntegerType& i) {
            if (i != 0)
                return (BIT_NUM - nlz64(i) - 1);
            else
                return 0; // arbitrary output for invalid argument 0
        }

        static IntegerType floorPow2(const IntegerType& i) {
            return flp2u64(i);
        }

        static IntegerType ceilPow2(const IntegerType& i) {
            return clp2u64(i);
        }
    };

    // ---------------------------------------------------------------
    // FUNCTIONS
    // ---------------------------------------------------------------

    template <typename I>
    I nlz(const I& i) {
        return MortonIntegerTraits<I>::nlz(i);
    }

    template <typename I>
    I log2Mod(const I& i) {
        return MortonIntegerTraits<I>::log2Mod(i);
    }

    /**
     * Given an input integer i, it returns:
     * 1) i > 0: the highest power of 2 less than or equal to i;
     * 2) i == 0: zero;
     * 3) i < 0: the opposite of the highest power of 2 less than or equal to i.
     * E.g.  floorPow2(6) = 4, floorPow2(0) = 0, floorPow2(-5) = -8.
     *        
     * @param i the input number
     * @return 
     */
    template <typename I>
    I floorPow2(const I& i) {
        return MortonIntegerTraits<I>::floorPow2(i);
    }

    /**
     * Given an input integer i, it returns:
     * 1) i > 0: the lowest power of 2 greater than or equal to i;
     * 2) i == 0: zero;
     * 3) i < 0: the opposite of the lowest power of 2 greater than or equal to i.
     * E.g.  ceilPow2(6) = 8, ceilPow2(0) = 0, ceilPow2(-5) = -4.
     *        
     * @param i the input number
     * @return 
     */
    template <typename I>
    I ceilPow2(const I& i) {
        return MortonIntegerTraits<I>::ceilPow2(i);
    }

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
    bool msb(const Integer& i1, const Integer& i2) {
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
    int mortonDistance(const Eigen::Matrix<I, Dim, 1>& v1, const Eigen::Matrix<I, Dim, 1>& v2) {
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

    template <typename I>
    int intervalPow2(const I& i1, const I& i2, I& low, I& mid, I& upp) {
        using MIT = MortonIntegerTraits<I>;
        using IntegerType = typename MIT::IntegerType;
        using UnsignedType = typename MIT::UnsignedType;
        IntegerType level, intervMask;

        level = MIT::BIT_NUM - MIT::nlz(MIT::removeSign(i1) ^ MIT::removeSign(i2));
        if (level < MIT::BIT_NUM) {
            intervMask = (1 << level) - 1;
            low = i1 & (~intervMask);
            upp = (low | intervMask);
            mid = low | (1 << (level - 1));
        } else {
            intervMask = ~0;
            low = std::numeric_limits<I>::min();
            upp = std::numeric_limits<I>::max();
            mid = (upp + low) / 2;
        }
        return (int)level;
    }

    template <typename I, int Dim>
    void mortonSplit(const Eigen::Matrix<I, Dim, 1>& v1, const Eigen::Matrix<I, Dim, 1>& v2,
            Eigen::Matrix<I, Dim, 1>& low, Eigen::Matrix<I, Dim, 1>& mid, Eigen::Matrix<I, Dim, 1>& upp) {
        int dimSplit, levelSplit, level;
        
        dimSplit = 0;
        levelSplit = intervalPow2(v1(0), v2(0), low(0), mid(0), upp(0));
//        ARS_VARIABLE2(dimSplit, levelSplit);
        for (int d = 1; d < Dim; ++d) {
            level = intervalPow2(v1(d), v2(d), low(d), mid(d), upp(d));
//            ARS_VARIABLE3(d, level, levelSplit);
            if (level > levelSplit) {
                mid(dimSplit) = low(dimSplit);
                dimSplit = d;
                levelSplit = level;
            }
            else {
                    mid(d) = low(d);
            }
        }
//        ARS_PRINT("split on dimension " << dimSplit << ", levelSplit " << levelSplit);
//        for (int d = 0; d < Dim; ++d) {
//            ARS_VARIABLE6(d, v1(d), v2(d), low(d), mid(d), upp(d));
//        }
    }

}

#endif /* MORTONSORT_H */

