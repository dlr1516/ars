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
#include <thrust/host_vector.h>
#include <map>
#include <algorithm>
#include <type_traits> // for static_assert()

#include <eigen3/Eigen/Dense>
#include <ars/definitions.h>

namespace ars {

    /**
     * Converts the input integer to its corresponding std::bitset. 
     * The function is designed for unsigned integer types, but also signed ones
     * can be used (after carefully thinking about the outcome).
     * @param i input integer
     * @param bs the std::bitset
     */
    template <typename Integer, unsigned int N>
    void convertIntegerToBitset(Integer i, std::bitset<N>& bs) {
        unsigned int bitNum = 8 * sizeof (Integer);
        bs.reset();
        for (unsigned int b = 0; b < N && b < bitNum; ++b) {
            bs[b] = i & (1 << b);
        }
    }

    /**
     * 
     * @param bsSrc
     * @param bsDst
     */
    template <unsigned int N, unsigned int M>
    void copyBitset(const std::bitset<N>& bsSrc, std::bitset<M>& bsDst) {
        bsDst.reset();
        for (unsigned int b = 0; b < N && b < M; ++b) {
            bsDst[b] = bsSrc[b];
        }
    }

    /**
     * Compares two std::bitsets with equal size when the two bitsets are 
     * interpreted as integer (bit 0 is the least significant one). 
     * Return true when bs1 is less than bs2. 
     * @param bs1 the first input bitset
     * @param bs2 the second input bitset
     */
    template <unsigned int N>
    bool lessBitsets(const std::bitset<N>& bs1, const std::bitset<N>& bs2) {
        for (int b = (int) N - 1; b >= 0; --b) {
            //ARS_VARIABLE4(b, bs1[b], bs2[b], bs1[b] ^ bs2[b]);
            if (bs1[b] ^ bs2[b]) return bs2[b];
        }
        return false;
    }

    template <unsigned int N>
    struct LessBitsets {

        bool operator()(const std::bitset<N>& bs1, const std::bitset<N>& bs2) const {
            return lessBitsets<N>(bs1, bs2);
        }
    };

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


                static const unsigned int SIMPLE_INDEX_BITNUM = Height; // number of bits to represent one coordinate
        static const unsigned int MULTI_INDEX_BITNUM = Dim * Height; // number of bits to represent Morton code
        static const unsigned int BRANCH_NUM = (1 << Dim); // number of children/branches of each node/octant

        using SimpleIndex = std::bitset<SIMPLE_INDEX_BITNUM>; // bitset for a single coordinate
        using MultiIndex = std::bitset<MULTI_INDEX_BITNUM>; // bitset for Morton code
        using ArraySimpleIndex = std::array<SimpleIndex, Dim>; // Dim-array of SimpleIndex
        using ArrayMultiIndex = std::array<MultiIndex, Dim>; // Dim-array of MultiIndex

        using Point = Eigen::Matrix<Scalar, Dim, 1>; // vector type to represent a point
        using VectorPoint = std::vector<Point, Eigen::aligned_allocator<Point> >;
        using MapPoint = std::multimap<
                MultiIndex,
                Point,
                LessBitsets<MULTI_INDEX_BITNUM>,
                Eigen::aligned_allocator<std::pair<MultiIndex, Point> >
                >;
        using Iterator = typename MapPoint::iterator;
        using ConstIterator = typename MapPoint::const_iterator;

        /**
         * Defaul constructor.
         */
        MortonOrderedPoints();

        /**
         * Destructor. 
         */
        ~MortonOrderedPoints();

        /**
         * Inserts the given vector of points and, if initBounds is set true, 
         * computes the minimum hypercube containing the input points. 
         * @param points the input point set
         * @param initBounds the flag that force computation of point region bounds
         */
        void insert(const VectorPoint& points, bool initBounds = true);

        /**
         * Returns the number of points in the data structure.
         * @return 
         */
        size_t size() const;

        /**
         * Returns the iterator to the first point according to Morton order.
         * @return 
         */
        ConstIterator begin() const;

        /**
         * Returns the iterator to the end of point set according to Morton order.
         * @return 
         */
        ConstIterator end() const;

        /**
         * Sets the bounds of the hypercube region of the point set
         * @param pmin vector of minimum coordinates
         * @param length size of the hypercube size
         */
        void setBounds(const Point& pmin, const Scalar& length);

        /**
         * Returns the maximum number of levels of the octree implicitly associated to 
         * MortonOrderedPoints. 
         */
        unsigned int getLevelMax() const;

        /**
         * Returns the number of octants for a given level. At a given level
         * each dimension is divided in 2^level bins and there are dim dimentions. 
         * Hence, octantNum = (2^level)^dim = 2^(level * dim)
         * E.g. l=1, d=2 (level 1 of quadtree) there are 4 octants (aka quadrants in two dimension) 
         *      l=2, d=2 (level 1 of quadtree) there are 16 octants, etc.
         * @param level
         * @return 
         */
        unsigned int getOctantNum(unsigned int level) const;

        /**
         * Returns the indices of the children nodes/octants of a given node/octant. 
         * @param octantParent
         * @param octantChildren
         */
        void getOctantChildren(unsigned int octantParent, std::array<unsigned int, BRANCH_NUM>& octantChildren) const;

        /**
         * Returns the minimum and maximum coordinates of the hyper-cube region 
         * delimiting an octant. 
         * @param level the level of the octree
         * @param octant the octant id 
         * @param pmin vector of minimum coordinates
         * @param pmax vector of maximum coordinates
         */
        void getOctantBounds(unsigned int level, unsigned int octant, Point& pmin, Point& pmax) const;

        /**
         * Returns the begin and end iterators to the points in the given hyper-cube 
         * octant. 
         * @param level the level of the octree
         * @param octant the octant id 
         * @param beg begin iterator
         * @param end end iterator
         */
        void getOctantPoints(unsigned int level, unsigned int octant, ConstIterator& beg, ConstIterator& end) const;

        /**
         * Encodes the input coordinate indices to an interlaced Morton code. 
         * @param indices the input array of indices (in bitset format)
         * @return 
         */
        static MultiIndex encode(const ArraySimpleIndex& indices);

        /**
         * Decode the interlaced Morton into an array of coordinate indices
         * @param mindex the input Morton interlaced code
         * @return 
         */
        static ArraySimpleIndex decode(const MultiIndex& mindex);

        //static MultiIndex enlarge(const SimpleIndex& si);

        //static SimpleIndex truncate(const MultiIndex& mi);

        /**
         * Converts a point in scalar coordinates to interlaced Morton code.
         * The conversion is based on discretization on Hypercube (supposedly) 
         * containing all the points 
         * @param p the point vector
         * @return 
         */
        inline MultiIndex pointToMorton(const Point& p) const;

        //        ConstIterator findLower(ConstIterator& beg, ConstIterator& end, const MultiIndex& miTarget) const;

        //        static bool compare(const MultiIndex& mi1, const MultiIndex& mi2);

    protected:
        //VectorPoint points_;
        MapPoint points_;
        Point pmin_;
        Scalar length_;
        uint64_t binNum_;
    };

    // ------------------------------------------------------------------------
    // ALIASES FOR MORTON TEMPLATE CLASSES
    // ------------------------------------------------------------------------

    using MortonOrderedPoints2f4h = MortonOrderedPoints<2, 4, float>;
    using MortonOrderedPoints2d4h = MortonOrderedPoints<2, 4, double>;
    using MortonOrderedPoints2f8h = MortonOrderedPoints<2, 8, float>;
    using MortonOrderedPoints2d8h = MortonOrderedPoints<2, 8, double>;
    using MortonOrderedPoints2f16h = MortonOrderedPoints<2, 16, float>;
    using MortonOrderedPoints2d16h = MortonOrderedPoints<2, 16, double>;
    using MortonOrderedPoints2f24h = MortonOrderedPoints<2, 24, float>;
    using MortonOrderedPoints2d24h = MortonOrderedPoints<2, 24, double>;
    using MortonOrderedPoints2f32h = MortonOrderedPoints<2, 32, float>;
    using MortonOrderedPoints2d32h = MortonOrderedPoints<2, 32, double>;

    using MortonOrderedPoints3f4h = MortonOrderedPoints<3, 4, float>;
    using MortonOrderedPoints3d4h = MortonOrderedPoints<3, 4, double>;
    using MortonOrderedPoints3f8h = MortonOrderedPoints<3, 8, float>;
    using MortonOrderedPoints3d8h = MortonOrderedPoints<3, 8, double>;
    using MortonOrderedPoints3f16h = MortonOrderedPoints<3, 16, float>;
    using MortonOrderedPoints3d16h = MortonOrderedPoints<3, 16, double>;
    using MortonOrderedPoints3f24h = MortonOrderedPoints<3, 24, float>;
    using MortonOrderedPoints3d24h = MortonOrderedPoints<3, 24, double>;
    using MortonOrderedPoints3f32h = MortonOrderedPoints<3, 32, float>;
    using MortonOrderedPoints3d32h = MortonOrderedPoints<3, 32, double>;

    // ------------------------------------------------------------------------
    // IMPLEMENTATION OF PUBLIC MEMBERS 
    // ------------------------------------------------------------------------

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    MortonOrderedPoints<Dim, Height, Scalar>::MortonOrderedPoints() : points_(), pmin_(Point::Zero()), length_(0) {
        static_assert(Dim > 0, "Dimension must be greater than or equal to 1");

        binNum_ = (uint64_t) (1) << (uint32_t) Height;
    }

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    MortonOrderedPoints<Dim, Height, Scalar>::~MortonOrderedPoints() {
    }

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    void MortonOrderedPoints<Dim, Height, Scalar>::insert(const VectorPoint& points, bool initBounds) {
        Point pmax;
        Scalar l;

        // Checks if the input point set is empty to avoid invalid access to memory 
        // later (we assume hence after that the point set has at least one point)
        if (points.empty()) {
            return;
        }

        // Copies the points into an internal vector and computes the min and max 
        // coordinates into pmin_ and pmax vectors
        //points_ = points;

        if (initBounds) {
            pmin_ = points[0];
            pmax = points[0];
            for (int i = 0; i < points.size(); ++i) {
                for (unsigned int d = 0; d < Dim; ++d) {
                    if (points[i](d) < pmin_(d))
                        pmin_(d) = points[i](d);
                    if (points[i](d) > pmax(d))
                        pmax(d) = points[i](d);
                }
            }
            ARS_VARIABLE2(pmin_.transpose(), pmax.transpose());

            // Computes the largest dimension
            length_ = 0.0;
            for (int d = 0; d < Dim; ++d) {
                l = pmax(d) - pmin_(d);
                if (l > length_)
                    length_ = l;
            }
            ARS_VARIABLE(length_);
        }

        for (auto it = points.begin(); it != points.end(); ++it) {
            std::cout << "  inserting point [" << it->transpose() << "] with code " << pointToMorton(*it)
                    << " points_.size() " << points_.size() << "\n";
            points_.insert(std::make_pair(pointToMorton(*it), *it));
        }

        ARS_PRINT("visiting the " << size() << " points inserted in Morton order:");
        for (auto it = points_.begin(); it != points_.end(); ++it) {
            std::cout << "  " << it->first << ": [" << it->second.transpose() << "]\n";
        }

        //        ARS_PRINT("point morton hash:");
        //        for (auto& p : points) {
        //            std::cout << "  " << p.transpose() << " -> " << pointToMorton(p) << "\n";
        //        }
        //
        // Sorts the points according to Morton order
        //        std::sort(points_.begin(), points_.end(),
        //                [&](const Point& p1, const Point & p2) -> bool {
        //                    //return pointToMorton(p1) < pointToMorton(p2);
        //                    //                    ARS_PRINT("compare p1 " << p1.transpose() << " (" << pointToMorton(p1) << ") "
        //                    //                            << ", p2 " << p2.transpose() << " (" << pointToMorton(p2) << "): "
        //                    //                            << "less(p1, p2) " << lessBitsets<MULTI_INDEX_BITNUM>(pointToMorton(p1), pointToMorton(p2)) << ", "
        //                    //                            << "less(p2, p1) " << lessBitsets<MULTI_INDEX_BITNUM>(pointToMorton(p2), pointToMorton(p1))
        //                    //                            );
        //                    //                    ARS_ASSERT(!(lessBitsets<MULTI_INDEX_BITNUM>(pointToMorton(p1), pointToMorton(p2)) && 
        //                    //                               lessBitsets<MULTI_INDEX_BITNUM>(pointToMorton(p2), pointToMorton(p1))));
        //                    return lessBitsets<MULTI_INDEX_BITNUM>(pointToMorton(p1), pointToMorton(p2));
        //                });
        //        ARS_PRINT("sorted pointset");
    }

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    void MortonOrderedPoints<Dim, Height, Scalar>::setBounds(const Point& pmin, const Scalar& length) {
        pmin_ = pmin;
        length_ = length;
    }

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    size_t
    MortonOrderedPoints<Dim, Height, Scalar>::size() const {
        return points_.size();
    }

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    typename MortonOrderedPoints<Dim, Height, Scalar>::ConstIterator
    MortonOrderedPoints<Dim, Height, Scalar>::begin() const {
        return points_.begin();
    }

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    typename MortonOrderedPoints<Dim, Height, Scalar>::ConstIterator
    MortonOrderedPoints<Dim, Height, Scalar>::end() const {
        return points_.end();
    }

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    unsigned int
    MortonOrderedPoints<Dim, Height, Scalar>::getLevelMax() const {
        return SIMPLE_INDEX_BITNUM;
    }

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    unsigned int
    MortonOrderedPoints<Dim, Height, Scalar>::getOctantNum(unsigned int level) const {
        //        unsigned int childrenNum = 2 << (Dim-1);
        //        unsigned int octantNum = 1;
        //
        //        for (int l = 0; l < level; ++l) {
        //            octantNum *= childrenNum;
        //        }
        //        //ARS_VARIABLE3(childrenNum, octantNum, level);
        //        return octantNum;
        return (1 << (Dim * level));
    }

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    void MortonOrderedPoints<Dim, Height, Scalar>::getOctantChildren(unsigned int octantParent, std::array<unsigned int, BRANCH_NUM>& octantChildren) const {
        octantChildren[0] = octantParent << Dim;
        for (unsigned int c = 1; c < BRANCH_NUM; ++c) {
            octantChildren[c] = octantChildren[0] + c;
        }
    }

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    void MortonOrderedPoints<Dim, Height, Scalar>::getOctantBounds(unsigned int level, unsigned int octant, Point& pmin, Point& pmax) const {
        ArraySimpleIndex indices;
        MultiIndex mi;
        Scalar octantLen;

        mi.reset();
        convertIntegerToBitset<unsigned int, MULTI_INDEX_BITNUM>(octant, mi);
        //mi = mi << (MULTI_INDEX_BITNUM - Dim * level);
        indices = decode(mi);

        octantLen = length_ / (1 << level);

        //        std::cout << "\n";
        //        ARS_VARIABLE2(mi, octantLen);
        for (int d = 0; d < Dim; ++d) {
            pmin(d) = pmin_(d) + octantLen * indices[d].to_ulong();
            pmax(d) = pmin(d) + octantLen;
            //            ARS_VARIABLE3(d, indices[d].to_ulong(), center(d));
        }

    }

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    void MortonOrderedPoints<Dim, Height, Scalar>::getOctantPoints(unsigned int level, unsigned int octant, ConstIterator& octBeg, ConstIterator& octEnd) const {
        MultiIndex miBeg, miEnd, miIt;
        unsigned int octantNum, octantBitStart;
        int count, step;
        ConstIterator octIt;
        //Point pOctMin, pOctMax;

        octantNum = getOctantNum(level);

        if (octant < 0 || octant >= octantNum) {
            std::cerr << __FILE__ << "," << __LINE__ << ": octant index " << octant << " overflow: "
                    << "on level " << level << " there are " << octantNum << " octants" << std::endl;
            return;
        }

        //        getOctantBounds(level, octant, pOctMin, pOctMax);
        //        beg = std::lower_bound(points_.begin(), points_.end(), pOctMin,
        //                [&](const Point& p1, const Point & p2) -> bool {
        //
        //                    return lessBitsets<MULTI_INDEX_BITNUM>(pointToMorton(p1), pointToMorton(p2));
        //                });
        //        end = std::lower_bound(points_.begin(), points_.end(), pOctMax,
        //                [&](const Point& p1, const Point & p2) -> bool {
        //
        //                    return lessBitsets<MULTI_INDEX_BITNUM>(pointToMorton(p1), pointToMorton(p2));
        //                });

        miBeg.reset();
        convertIntegerToBitset<unsigned int, MULTI_INDEX_BITNUM>(octant, miBeg);
        miBeg = miBeg << (MULTI_INDEX_BITNUM - Dim * level);

        miEnd.set();
        miEnd = (miEnd >> (Dim * level)) | miBeg;

        octBeg = points_.lower_bound(miBeg);
        octEnd = points_.upper_bound(miEnd);

        //ARS_VARIABLE4(miBeg, octBeg->first, miEnd, octEnd->first);

        //octBeg = findLower(points_.begin(), points_.end(), miBeg);
        //octEnd = findLower(octBeg, std::end(points_), miEnd);

        //        ARS_VARIABLE4(miBeg, pointToMorton(pOctMin), miEnd, pointToMorton(pOctMax));
        //
        //        beg = std::find_if(points_.begin(), points_.end(),
        //                [&](const Point & p) {
        //                    return !lessBitsets<MULTI_INDEX_BITNUM>(pointToMorton(p), miBeg); });
        //        end = std::find_if(points_.begin(), points_.end(),
        //                [&](const Point & p) {
        //                    return lessBitsets<MULTI_INDEX_BITNUM>(pointToMorton(p), miEnd); });
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

    //    template <unsigned int Dim, unsigned int Height, typename Scalar>
    //    typename MortonOrderedPoints<Dim, Height, Scalar>::MultiIndex
    //    MortonOrderedPoints<Dim, Height, Scalar>::enlarge(const SimpleIndex& si) {
    //        MultiIndex mi;
    //        mi.reset();
    //        for (unsigned int b = 0; b < SIMPLE_INDEX_BITNUM; ++b) {
    //            mi[b] = si[b];
    //        }
    //        return mi;
    //    }

    //    template <unsigned int Dim, unsigned int Height, typename Scalar>
    //    typename MortonOrderedPoints<Dim, Height, Scalar>::SimpleIndex
    //    MortonOrderedPoints<Dim, Height, Scalar>::truncate(const MultiIndex& mi) {
    //        SimpleIndex si;
    //        for (unsigned int b = 0; b < SIMPLE_INDEX_BITNUM; ++b) {
    //            si[b] = mi[b];
    //        }
    //        return si;
    //    }

    template <unsigned int Dim, unsigned int Height, typename Scalar>
    typename MortonOrderedPoints<Dim, Height, Scalar>::MultiIndex
    MortonOrderedPoints<Dim, Height, Scalar>::pointToMorton(const Point& p) const {
        ArraySimpleIndex indices;
        unsigned int coordIndex;

        for (int d = 0; d < Dim; ++d) {
            coordIndex = (unsigned int) floor((p(d) - pmin_(d)) / length_ * binNum_);
            if (coordIndex >= binNum_) {
                coordIndex = binNum_ - 1;
            }
            convertIntegerToBitset<unsigned int, SIMPLE_INDEX_BITNUM>(coordIndex, indices[d]);
            //ARS_VARIABLE4(d, p(d), coordIndex, indices[d]);
        }
        //ARS_VARIABLE(encode(indices));
        return encode(indices);
    }

    //    template <unsigned int Dim, unsigned int Height, typename Scalar>
    //    typename MortonOrderedPoints<Dim, Height, Scalar>::ConstIterator
    //    MortonOrderedPoints<Dim, Height, Scalar>::findLower(ConstIterator& beg, ConstIterator& end, const MultiIndex& miTarget) const {
    //        ConstIterator it;
    //        MultiIndex miIt;
    //        typename std::iterator_traits<ConstIterator>::difference_type count, step;
    //        count = std::distance(beg, end);
    //
    //        while (count > 0) {
    //            it = beg;
    //            step = count / 2;
    //            std::advance(it, step);
    //            miIt = pointToMorton(*it);
    //            if (!lessBitsets<MULTI_INDEX_BITNUM>(miTarget, miIt)) {
    //                beg = ++it;
    //                count -= step + 1;
    //            } else
    //                count = step;
    //        }
    //        return beg;
    //    }

    //    template <unsigned int Dim, unsigned int Height, typename Scalar>
    //    bool MortonOrderedPoints<Dim, Height, Scalar>::compare(const MultiIndex& mi1, const MultiIndex& mi2) {
    //        for (int b = MULTI_INDEX_BITNUM - 1; b >= 0; --b) {
    //            if (mi1[b] ^ mi2[b]) return mi2[b];
    //        }
    //        return false;
    //    }

} // end of namespace

#endif /* MORTONORDEREDPOINTS_H */

