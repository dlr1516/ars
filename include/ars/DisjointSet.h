/**
 * ARS - Angular Radon Spectrum 
 * Copyright (C) 2017 Dario Lodi Rizzini.
 *           (C) 2021 Dario Lodi Rizzini, Ernesto Fontana.
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
#ifndef DISJOINT_SET_H
#define DISJOINT_SET_H

#include <iostream>
#include <vector>

namespace ars {

    /**
     * Class DisjointSet is a lightweight implementation of the standard Robert Tarjan 
     * Disjoint-Set data structure.
     * 
     * Tarjan, R.E. "Efficiency of a Good But Not Linear Set Union Algorithm". 
     * Journal of the ACM 22 (2), 215–225, 1975. doi:10.1145/321879.321884.
     */
    class DisjointSet {
    public:
        typedef int id_type;

        struct Node {
            id_type parent;
            id_type rank;
            size_t size;
        };

        /** Creates an empty disjoint set.
         */
        DisjointSet();

        /** Creates a disjoint set with an initial number of nodes.
         * Each node corresponds to a set at the beginning. 
         */
        DisjointSet(id_type num);

        /**
         * Destructor. 
         */
        ~DisjointSet();

        /** Creates a disjoint set with an initial number of nodes.
         * Each node corresponds to a set at the beginning. 
         */
        void initSets(id_type num);

        /** Creates a new node corresponding to a set (with one item). 
         */
        id_type makeSet();

        /** Returns the equivalent set identifier. 
         */
        id_type find(id_type n);

        /** Joins two sets.
         */
        id_type join(id_type n1, id_type n2);

        /**
         * Returns the number of nodes.
         */
        size_t nodeNum() const {
            return nodes_.size();
        }

        /**
         * Returns the number of nodes.
         */
        size_t size() const {
            return nodeNum();
        }

        /** Returns the size of cluster with the given label.
         */
        size_t size(int n) {
            id_type root = find(n);
            return nodes_[root].size;
        }

        /** Returns the number of label sets. 
         */
        int setNum() const {
            return setNum_;
        }

        /** Returns the list of parents.
         */
        template <typename InsertIt>
        void parents(InsertIt it) {
            std::vector<Node>::iterator nit;
            id_type id = 0;
            for (nit = nodes_.begin(); nit != nodes_.end(); ++nit, ++id) {
                if (nit->parent == id) {
                    it = id;
                }
            }
        }

    private:
        std::vector<Node> nodes_;
        int setNum_;
    };

}

#endif

