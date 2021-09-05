#ifndef MORTONOCTREE_H
#define MORTONOCTREE_H

#include <iostream>
#include <cstdint> 
#include <vector>

#include <Eigen/Dense>
#include <ars/definitions.h>
#include <ars/MortonSort.h>

namespace ars {

    template <unsigned int Dim, typename Scalar, typename Integer = int64_t>
    class MortonOctree {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        using PointInt = Eigen::Matrix<Integer, Dim, 1>;
        using PointScalar = Eigen::Matrix<Scalar, Dim, 1>;
        using VectorPointInt = std::vector<PointInt, Eigen::aligned_allocator<PointInt> >;
        using VectorPointScalar = std::vector<PointScalar, Eigen::aligned_allocator<PointScalar> >;

        struct Item {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            PointInt index;
            PointScalar value;
            size_t pos;
        };
        using VectorItem = std::vector<Item, Eigen::aligned_allocator<Item> >;
        using Iterator = typename VectorItem::iterator;
        using ConstIterator = typename VectorItem::const_iterator;

        MortonOctree() : items_(), res_(1.0) {
        }

        virtual ~MortonOctree() {
        }

        void setRes(Scalar res) {
            res_ = res;
        }

        void insert(const VectorPointScalar& points) {
            // Converts each input point into an item (original point, integer 
            // indices, position in original vector)
            items_.resize(points.size());
            for (size_t i = 0; i < points.size(); ++i) {
                pointToItem(points[i], i, items_[i]);
            }

            // Sorts in morton order
            std::sort(std::begin(items_), std::end(items_),
                    [&](const Item& item1, const Item & item2) -> bool {
                        return mortonCmpInt(item1.index, item2.index);
                    }
            );
        }

        const VectorItem& getItems() const {
            return items_;
        }

        size_t size() const {
            return items_.size();
        }

        Iterator begin() {
            return items_.begin();
        }

        Iterator end() {
            return items_.end();
        }

        ConstIterator begin() const {
            return items_.begin();
        }

        ConstIterator end() const {
            return items_.end();
        }

        ConstIterator find(const PointScalar& q) const {
            Item item;
            pointToItem(q, 0, item);
            ConstIterator it = std::lower_bound(std::begin(items_), std::end(items_), item,
                    [&](const Item& item1, const Item & item2) -> bool {
                        return mortonCmpInt(item1.index, item2.index);
                    }
            );
            return it;
        }


    private:
        VectorItem items_;
        Scalar res_;

        void pointToItem(const PointScalar& p, size_t pos, Item& item) const {
            item.value = p;
            for (size_t d = 0; d < Dim; ++d) {
                item.index(d) = (Integer) round(p(d) / res_);
            }
            item.pos = pos;
        }
    };

}

#endif /* MORTONOCTREE_H */

