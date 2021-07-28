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
#ifndef MORTON_TREE_H
#define MORTON_TREE_H

#include <iostream>
#include <cstdint> 
#include <array>
#include <vector>
#include <map>
#include <algorithm>
#include <type_traits> // for static_assert()

#include <Eigen/Dense>
#include <ars/definitions.h>
#include <bitset>

namespace ars {

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

    template <unsigned int Dim, typename Item, typename Integer = uint64_t>
    class MortonTree {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        static const Integer INTEGER_BITNUM = std::numeric_limits<Integer>::digits;
        static const Integer INTEGER_MAX = std::numeric_limits<Integer>::max();
        static const Integer CHILDREN_NUM = 1 << Dim;

        using IntPoint = Eigen::Matrix<Integer, Dim, 1>;

        struct Octant {
            int level;
            IntPoint mask;

            std::string getString() const {
                std::stringstream os;
                os << "(level " << level << " [";
                for (int d = 0; d < Dim; ++d) {
                    os << (mask(d) >> (INTEGER_BITNUM - level)) << " ";
                }
                os << "])";
                return os.str();
            }
        };
        using OctanChildrenArray = std::array<Octant, CHILDREN_NUM>;

        struct LessIntPoint {
            Integer shift;

            LessIntPoint();

            LessIntPoint(Integer s);

            bool operator()(const IntPoint& p1, const IntPoint& p2) const;
        };

        using MapPoint = std::multimap<
                IntPoint,
                Item,
                LessIntPoint,
                Eigen::aligned_allocator<std::pair<IntPoint, Item> >
                >;
        using Iterator = typename MapPoint::iterator;
        using ConstIterator = typename MapPoint::const_iterator;


        /**
         * Default constructor.
         */
        MortonTree();

        /**
         * Destructor. 
         */
        ~MortonTree();

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
         * Returns the maximum number of levels of the octree implicitly associated to 
         * MortonOrderedPoints. 
         */
        unsigned int getLevelMax() const;

        unsigned int getLevel(const Octant& octant) const;

        Octant root() const;

        Octant parent(const Octant& octant) const;

        OctanChildrenArray children(const Octant& octant) const;

        void getOctant(const Octant& octant, ConstIterator& octantBeg, ConstIterator& octantEnd) const;

        /**
         * Inserts a (integer) point into the Morton tree with its corresponding Item. 
         * @param point the point inserted
         * @param item
         */
        void insert(const IntPoint& point, const Item& item);


        /**
         */
        //void getOctant(const Octant& octant, ConstIterator& octantBeg, ConstIterator& octantEnd) const;



    private:
        MapPoint data_;
    };


    // ------------------------------------------------------------------------
    // IMPLEMENTATION OF LessIntPoint
    // ------------------------------------------------------------------------

    template <unsigned int Dim, typename Item, typename Integer>
    MortonTree<Dim, Item, Integer>::LessIntPoint::LessIntPoint() : shift(0) {
    }

    template <unsigned int Dim, typename Item, typename Integer>
    MortonTree<Dim, Item, Integer>::LessIntPoint::LessIntPoint(Integer s) : shift(s) {
    }

    template <unsigned int Dim, typename Item, typename Integer>
    bool MortonTree<Dim, Item, Integer>::LessIntPoint::operator()(const IntPoint& p1, const IntPoint& p2) const {
        Integer dimLast = 0;
        Integer diffLast = 0;
        Integer diffShifted;
        for (Integer d = 0; d < Dim; ++d) {
            diffShifted = (p1(d) + shift) ^ (p2(d) + shift);
            if (lessMbs(diffLast, diffShifted)) {
                dimLast = d;
                diffLast = diffShifted;
            }
        }
        return p1(dimLast) < p2(dimLast);
    }

    // ------------------------------------------------------------------------
    // IMPLEMENTATION OF MortonTree
    // ------------------------------------------------------------------------

    template <unsigned int Dim, typename Item, typename Integer>
    MortonTree<Dim, Item, Integer>::MortonTree() : data_() {
    }

    template <unsigned int Dim, typename Item, typename Integer>
    MortonTree<Dim, Item, Integer>::~MortonTree() {
    }

    template <unsigned int Dim, typename Item, typename Integer>
    size_t MortonTree<Dim, Item, Integer>::size() const {
        return data_.size();
    }

    template <unsigned int Dim, typename Item, typename Integer>
    typename MortonTree<Dim, Item, Integer>::ConstIterator
    MortonTree<Dim, Item, Integer>::begin() const {
        return data_.begin();
    }

    template <unsigned int Dim, typename Item, typename Integer>
    typename MortonTree<Dim, Item, Integer>::ConstIterator
    MortonTree<Dim, Item, Integer>::end() const {
        return data_.end();
    }

    template <unsigned int Dim, typename Item, typename Integer>
    unsigned int MortonTree<Dim, Item, Integer>::getLevelMax() const {
        return INTEGER_BITNUM;
    }

    template <unsigned int Dim, typename Item, typename Integer>
    typename MortonTree<Dim, Item, Integer>::Octant
    MortonTree<Dim, Item, Integer>::root() const {
        Octant oct;
        oct.level = 0;
        oct.mask = IntPoint::Zero();
        return oct;
    }

    template <unsigned int Dim, typename Item, typename Integer>
    typename MortonTree<Dim, Item, Integer>::Octant
    MortonTree<Dim, Item, Integer>::parent(const Octant& octant) const {
        Octant par;
        par.level = octant.level - 1;
        for (int d = 0; d < Dim; ++d) {
            par.mask(d) = (octant.mask(d) >> 1);
        }
        return par;
    }

    template <unsigned int Dim, typename Item, typename Integer>
    typename MortonTree<Dim, Item, Integer>::OctanChildrenArray
    MortonTree<Dim, Item, Integer>::children(const Octant& octant) const {
        std::array<Octant, CHILDREN_NUM> children;
        Integer digit;
        for (int c = 0; c < CHILDREN_NUM; ++c) {
            children[c].level = octant.level + 1;
            //            ARS_PRINT("child node: c " << c << " level" << children[c].level);
            for (int d = 0; d < Dim; ++d) {
                digit = c >> d;
                children[c].mask(d) = octant.mask(d) | (digit << (INTEGER_BITNUM - children[c].level));
                //                std::cout << "  d " << d << ": children[c].mask(d) " << children[c].mask(d) 
                //                        << " (" << std::bitset<INTEGER_BITNUM>(children[c].mask(d)) << ")\n";
            }
        }
        return children;
    }

    template <unsigned int Dim, typename Item, typename Integer>
    void MortonTree<Dim, Item, Integer>::getOctant(const Octant& octant, ConstIterator& octantBeg, ConstIterator& octantEnd) const {
        IntPoint pmin, pmax;
        Integer mask = (INTEGER_MAX >> octant.level);

        ARS_PRINT("bounds octant: level " << octant.level << " mask " << octant.mask);
        pmin = octant.mask;
        for (int d = 0; d < Dim; ++d) {
            pmax(d) = pmin(d) | mask;
            std::cout << "  d " << d << ": pmin " << pmin(d) << " (" << std::bitset<INTEGER_BITNUM>(pmin(d)) << ") "
                    << "pmax " << pmax(d) << " (" << std::bitset<INTEGER_BITNUM>(pmax(d)) << ")\n";
        }

        octantBeg = data_.lower_bound(pmin);
        octantEnd = data_.upper_bound(pmax);
                
        ARS_PRINT("list of octant items");
        for (auto it = octantBeg; it != octantEnd; ++it) {
            ARS_VARIABLE2(it->first, it->second);
        }
    }

    template <unsigned int Dim, typename Item, typename Integer>
    void MortonTree<Dim, Item, Integer>::insert(const IntPoint& point, const Item& item) {
        data_.insert(std::make_pair(point, item));
    }

    // ------------------------------------------------------------------------
    // OTHERS
    // ------------------------------------------------------------------------

    //    template <unsigned int Dim, typename Item, typename Integer>
    //    std::ostream& operator<<(std::ostream& os, const typename MortonTree<Dim, Item, Integer>::Octant& octant) {
    //        os << "(level " << octant.level << " [";
    //        for (int d = 0; d < Dim; ++d) {
    //            os << (octant.mask >> (MortonTree<Dim, Item, Integer>::INTEGER_BITNUM - octant.level)) << " ";
    //        }
    //         os << "])";
    //         return os;
    //     }
};

#endif /* MORTON_TREE_H */

