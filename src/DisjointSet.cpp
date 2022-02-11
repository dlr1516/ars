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
#include <ars/DisjointSet.h>

namespace cuars {

    DisjointSet::DisjointSet() : setNum_(0) {
    }

    DisjointSet::DisjointSet(id_type num) : setNum_(num) {
        Node n;
        for (id_type i = 0; i < num; ++i) {
            n.parent = i;
            n.rank = 0;
            n.size = 1;
            nodes_.push_back(n);
        }
    }

    DisjointSet::~DisjointSet() {
    }

    void DisjointSet::initSets(id_type num) {
        Node n;
        nodes_.clear();
        setNum_ = num;
        for (id_type i = 0; i < num; ++i) {
            n.parent = i;
            n.rank = 0;
            n.size = 1;
            nodes_.push_back(n);
        }
    }

    DisjointSet::id_type DisjointSet::makeSet() {
        Node n;
        n.parent = (id_type) nodes_.size();
        n.rank = 0;
        n.size = 1;
        nodes_.push_back(n);
        return n.parent;
    }

    DisjointSet::id_type DisjointSet::find(id_type n) {
        if (n < 0 || n >= (id_type) nodes_.size()) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid id " << n
                    << " should be positive less than " << nodes_.size() << "!" << std::endl;
            return nodes_.size();
        }
        id_type tmp = n;
        while (tmp != nodes_[tmp].parent) {
            tmp = nodes_[tmp].parent;
        }
        nodes_[n].parent = tmp; // path compression: reduces the hierarchical rank
        return tmp;
    }

    /** Joins two sets.
     */
    DisjointSet::id_type DisjointSet::join(id_type n1, id_type n2) {
        id_type root1 = find(n1);
        id_type root2 = find(n2);
        // If n1 and n2 are already in the same set, there is no join to do
        if (root1 == root2) return nodes_.size();

        setNum_--;
        if (nodes_[root1].rank < nodes_[root2].rank) {
            nodes_[root1].parent = root2;
            nodes_[root2].size += nodes_[root1].size;
            return root2;
        } else if (nodes_[root1].rank > nodes_[root2].rank) {
            nodes_[root2].parent = root1;
            nodes_[root1].size += nodes_[root2].size;
            return root1;
        } else {
            nodes_[root2].parent = root1;
            nodes_[root1].size += nodes_[root2].size;
            nodes_[root1].rank++;
            return root1;
        }
    }

}

