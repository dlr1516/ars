/**
 * ARS - Angular Radon Spectrum 
 * Copyright (C) 2017-2020 Dario Lodi Rizzini.
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
#ifndef MORTONORDEREDPOINTS_H
#define MORTONORDEREDPOINTS_H

#include <iostream>
#include <cstdint>
#include <bitset>
#include <array>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <ars/definitions.h>

namespace ars {
    
    template <typename Integer, unsigned int N>
    void convertIntegerToBitset(Integer i, std::bitset<N>& bs) {
        unsigned int bitNum = 8 * sizeof(Integer);
        bs.reset();
        for (unsigned int b = 0; b < N && b < bitNum; ++b) {
            bs[b] = i & (1 << b);
        }
    }
    
    template <unsigned int N>
    bool compareBitsets(const std::bitset<N>& bs1, const std::bitset<N>& bs2) {
        for (unsigned int b = N - 1; b >= 0; --b) {
            if (bs1[b] ^ bs2[b]) return bs2[b];
        }
        return false;
    }

    /**
     * Class MortonOrderedPoints provides a simple data structure for sorting the points 
     * of an input point set (fixed) according to Morton order. 
     * The parameter of the class are:
     * - Dim: the dimension of the space;
     * - Height: the maximum height of the octree associated to the Morton order indexing,
     *   i.e. the number of intervals each axis is divided or equivalently 
     *   the number of bits used to represent the discretized coordinate;
     * - Scalar: the scalar type used to represent the point coordinate. 
     */
    template <unsigned int Dim, unsigned int Height, typename Scalar = double>
    class MortonOrderedPoints {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

                static const unsigned int SIMPLE_INDEX_BITNUM = Height;
        static const unsigned int MULTI_INDEX_BITNUM = Dim * Height;

        using SimpleIndex = std::bitset<SIMPLE_INDEX_BITNUM>;
        using MultiIndex = std::bitset<MULTI_INDEX_BITNUM>;
        using ArraySimpleIndex = std::array<SimpleIndex, Dim>;
        using ArrayMultiIndex = std::array<MultiIndex, Dim>;

        using Point = Eigen::Matrix<Scalar, Dim, 1>;
        using VectorPoint = std::vector<Point, Eigen::aligned_allocator<Point> >;
        using Iterator = typename VectorPoint::iterator;
        using ConstIterator = typename VectorPoint::const_iterator;

        MortonOrderedPoints();

        ~MortonOrderedPoints();

        unsigned int getOctantNum(unsigned int level) const;

        void insert(const VectorPoint& points);

        void getInterval(unsigned int level, unsigned int quadrant, ConstIterator& beg, ConstIterator& end) const;

        static MultiIndex encode(const ArraySimpleIndex& indices);

        static ArraySimpleIndex decode(const MultiIndex& mindex);
        
        static MultiIndex enlarge(const SimpleIndex& si);

        static SimpleIndex truncate(const MultiIndex& mi);

        inline MultiIndex pointToMorton(const Point& p) const;

//        static bool compare(const MultiIndex& mi1, const MultiIndex& mi2);

    protected:
        VectorPoint points_;
        Point pmin_;
        Scalar length_;
    };
    
    using MortonOrderedPoints2f4h = MortonOrderedPoints<2, 4, float>;
    using MortonOrderedPoints2d4h = MortonOrderedPoints<2, 4, double>;
    using MortonOrderedPoints2f8h = MortonOrderedPoints<2, 8, float>;
    using MortonOrderedPoints2d8h = MortonOrderedPoints<2, 8, double>;
    using MortonOrderedPoints2f16h = MortonOrderedPoints<2, 16, float>;
    using MortonOrderedPoints2d16h = MortonOrderedPoints<2, 16, double>;
    using MortonOrderedPoints2f24h = MortonOrderedPoints<2, 24, float>;
    using MortonOrderedPoints2d24h = MortonOrderedPoints<2, 24, double>;

    // ------------------------------------------------------------------------
    // IMPLEMENTATION OF PUBLIC MEMBERS 
    // ------------------------------------------------------------------------

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    MortonOrderedPoints<Dim, Height, Scalar>::MortonOrderedPoints() : points_(), pmin_(Point::Zero()), length_(0) {
    }

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    MortonOrderedPoints<Dim, Height, Scalar>::~MortonOrderedPoints() {
    }

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    unsigned int
    MortonOrderedPoints<Dim, Height, Scalar>::getOctantNum(unsigned int level) const {
        unsigned int childrenNum = 2 << Dim;
        unsigned int octantNum = 1;

        for (int l = 0; l < level; ++l) {
            octantNum *= childrenNum;
        }
        return octantNum;
    }

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    void MortonOrderedPoints<Dim, Height, Scalar>::insert(const VectorPoint& points) {
        Point pmax;
        Scalar l;

        // Checks if the input point set is empty to avoid invalid access to memory 
        // later (we assume hence after that the point set has at least one point)
        if (points.empty()) {
            return;
        }

        // Copies the points into an internal vector and computes the min and max 
        // coordinates into pmin_ and pmax vectors
        points_ = points;
        pmin_ = points_[0];
        pmax = points_[0];
        for (int i = 0; i < points_.size(); ++i) {
            for (unsigned int d = 0; d < Dim; ++d) {
                if (points_[i](d) < pmin_(d))
                    pmin_[i] = points_[i](d);
                if (points_[i](d) > pmax(d))
                    pmax(d) = points_[i](d);
            }
        }

        // Computes the largest dimension
        length_ = 0.0;
        for (int d = 0; d < Dim; ++d) {
            l = pmax(d) - pmin_(d);
            if (l > length_)
                length_ = l;
        }

        // Sorts the points according to Morton order
        std::sort(points_.begin(), points_.end(),
                [&](const Point& p1, const Point & p2) -> bool {
                    //return pointToMorton(p1) < pointToMorton(p2);
                    return compareBitsets<MULTI_INDEX_BITNUM>(pointToMorton(p1), pointToMorton(p2));
                });
    }

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    void MortonOrderedPoints<Dim, Height, Scalar>::getInterval(unsigned int level, unsigned int octant, ConstIterator& beg, ConstIterator& end) const {
        MultiIndex miBeg, miEnd; //, prefix;
        unsigned int octantNum, octantBitStart;
        
        octantNum = getOctantNum(level);

        if (octant < 0 || octant >= octantNum) {
            std::cerr << __FILE__ << "," << __LINE__ << ": octant index " << octant << " overflow:"
                    << "on level " << level << " there are " << octantNum << " octants" << std::endl;
            return;
        }

        miBeg.reset();
        miEnd.set();
        //prefix.to_ulong() = octant;
        octantBitStart = MULTI_INDEX_BITNUM - 1 - level;
        for (unsigned int b = MULTI_INDEX_BITNUM - 1; b >= octantBitStart; --b) {
            if (b < level) {
//                miBeg[b] = prefix[MULTI_INDEX_BITNUM - 1 - b];
//                miEnd[b] = prefix[MULTI_INDEX_BITNUM - 1 - b];
                miBeg[b] = (octant & (1 << (b - octantBitStart)));
                miEnd[b] = miBeg[b];
            }
        }

        beg = std::find_if(points_.begin(), points_.end(),
                [&](const Point & p) {
                    return compareBitsets<MULTI_INDEX_BITNUM>(pointToMorton(p), miBeg); });
        end = std::find_if(points_.begin(), points_.end(),
                [&](const Point & p) {
                    return compareBitsets<MULTI_INDEX_BITNUM>(pointToMorton(p), miEnd); });
    }

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    typename MortonOrderedPoints<Dim, Height, Scalar>::MultiIndex
    MortonOrderedPoints<Dim, Height, Scalar>::encode(const ArraySimpleIndex& indices) {
        MultiIndex res;
//        MultiIndex one;

        res.reset();
//        one.reset();
//        one.set(0) = 1;
        for (unsigned int b = 0; b < SIMPLE_INDEX_BITNUM; ++b) {
            for (unsigned int d = 0; d < Dim; ++d) {
                //res |= (enlarge(indices[d]) & (one << b)) << b;
                res[Dim * b + d] = indices[d][b];
            }
        }

        return res;
    }

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    typename MortonOrderedPoints<Dim, Height, Scalar>::ArraySimpleIndex
    MortonOrderedPoints<Dim, Height, Scalar>::decode(const MultiIndex& mindex) {
        ArraySimpleIndex indices;
//        MultiIndex one;

        // Resets the indices
//        one.reset();
//        one.set(0) = 1;
//        for (unsigned int d = 0; d < Dim; ++d) {
//            indices[d].reset();
//        }

        for (unsigned int b = 0; b < SIMPLE_INDEX_BITNUM; ++b) {
            for (unsigned int d = 0; d < Dim; ++d) {
//                indices[d] |= truncate((mindex & (one << (Dim * b + d))) >> (b + d));
                indices[d][b] = mindex[Dim * b + d];
            }
        }
        return indices;
    }

    // ------------------------------------------------------------------------
    // IMPLEMENTATION OF PRIVATE/PROTECTED MEMBERS 
    // ------------------------------------------------------------------------

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    typename MortonOrderedPoints<Dim, Height, Scalar>::MultiIndex
    MortonOrderedPoints<Dim, Height, Scalar>::enlarge(const SimpleIndex& si) {
        MultiIndex mi;
        mi.reset();
        for (unsigned int b = 0; b < SIMPLE_INDEX_BITNUM; ++b) {
            mi[b] = si[b];
        }
        return mi;
    }

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    typename MortonOrderedPoints<Dim, Height, Scalar>::SimpleIndex
    MortonOrderedPoints<Dim, Height, Scalar>::truncate(const MultiIndex& mi) {
        SimpleIndex si;
        for (unsigned int b = 0; b < SIMPLE_INDEX_BITNUM; ++b) {
            si[b] = mi[b];
        }
        return si;
    }

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    typename MortonOrderedPoints<Dim, Height, Scalar>::MultiIndex
    MortonOrderedPoints<Dim, Height, Scalar>::pointToMorton(const Point& p) const {
        ArraySimpleIndex indices;
        unsigned int coordIndex;
        for (int d = 0; d < Dim; ++d) {
            coordIndex = (unsigned int) floor((p(d) - pmin_(d)) / length_ * SIMPLE_INDEX_BITNUM);
            convertIntegerToBitset<unsigned int, SIMPLE_INDEX_BITNUM>(coordIndex, indices[d]);
        }
        return encode(indices);
    }

//    template <unsigned int Dim, unsigned int Height, typename Scalar>
//    bool MortonOrderedPoints<Dim, Height, Scalar>::compare(const MultiIndex& mi1, const MultiIndex& mi2) {
//        for (int b = MULTI_INDEX_BITNUM - 1; b >= 0; --b) {
//            if (mi1[b] ^ mi2[b]) return mi2[b];
//        }
//        return false;
//    }

} // end of namespace

#endif /* MORTONORDEREDPOINTS_H */

