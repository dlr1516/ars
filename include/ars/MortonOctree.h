#ifndef MORTONOCTREE_H
#define MORTONOCTREE_H

#include <cstdint>
#include <iostream>
#include <vector>

#include <Eigen/Dense>
//#include <ars/definitions.h>
//#include <ars/MortonSort.h>
#include "MortonSort.h"
#include "definitions.h"

namespace ars {

/**
 * Class MortonOctree is a data structure that stores scalar points in
 * a vector according to Morton order.
 * Morton-ordered vectors implicitly represent quad/octrees.
 */
template <unsigned int Dim, typename Scalar, typename Integer = int64_t>
class MortonOctree {
   public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

#if __cplusplus < 201703L
    using PointInt = Eigen::Matrix<Integer, Dim, 1>;
    using PointScalar = Eigen::Matrix<Scalar, Dim, 1>;
    using VectorPointInt =
        std::vector<PointInt, Eigen::aligned_allocator<PointInt> >;
    using VectorPointScalar =
        std::vector<PointScalar, Eigen::aligned_allocator<PointScalar> >;

    struct Item {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        PointInt index;
        PointScalar value;
        size_t pos;
    };
    using VectorItem = std::vector<Item, Eigen::aligned_allocator<Item> >;
#else
    using PointInt = Eigen::Matrix<Integer, Dim, 1>;
    using PointScalar = Eigen::Matrix<Scalar, Dim, 1>;
    using VectorPointInt = std::vector<PointInt>;
    using VectorPointScalar = std::vector<PointScalar>;

    struct Item {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        PointInt index;
        PointScalar value;
        size_t pos;
    };
    using VectorItem = std::vector<Item>;
#endif

    using Iterator = typename VectorItem::iterator;
    using ConstIterator = typename VectorItem::const_iterator;

    /**
     * Default constructor.
     */
    MortonOctree() : items_(), res_(1.0) {}

    /**
     * Destructor.
     */
    virtual ~MortonOctree() {}

    /**
     * Sets the minimum resolution of octree.
     * The whole Dim-dimension space is dicretized into cells with the given
     * resolution.
     * @param res the value of resolution
     */
    void setRes(Scalar res) { res_ = res; }

    /**
     * Inserts the points of the input vector into the octree.
     * The points are internally stored inside struct Item and
     * sorted according to Morton order (aka z-order).
     * @param points the input points
     */
    void insert(const VectorPointScalar& points) {
        // Converts each input point into an item (original point, integer
        // indices, position in original vector)
        items_.resize(points.size());
        for (size_t i = 0; i < points.size(); ++i) {
            pointToItem(points[i], i, items_[i]);
        }

        // Sorts in morton order
        std::sort(std::begin(items_), std::end(items_),
                  [&](const Item& item1, const Item& item2) -> bool {
                      return mortonCmpInt(item1.index, item2.index);
                  });
    }

    /**
     * Returns a reference to the internal container.
     * @return
     */
    const VectorItem& getItems() const { return items_; }

    /**
     * Returns the number of items.
     * @return
     */
    size_t size() const { return items_.size(); }

    /**
     * Returns the iterator pointing to the first item of the octree vector.
     * @return
     */
    Iterator begin() { return items_.begin(); }

    /**
     * Returns the iterator pointing to the ending position of the octree
     * vector. (Ending position is after the last item).
     * @return
     */
    Iterator end() { return items_.end(); }

    /**
     * Returns the iterator pointing to the first item of the octree vector.
     * @return
     */
    ConstIterator begin() const { return items_.begin(); }

    /**
     * Returns the iterator pointing to the ending position of the octree
     * vector. (Ending position is after the last item).
     * @return
     */
    ConstIterator end() const { return items_.end(); }

    /**
     * Returns the iterator pointing to the item nearest to the query point.
     * The nearest point is the item NOT LESS than the query item in Morton
     * order.
     * @param q the query point
     * @return
     */
    ConstIterator findNearest(const PointScalar& q) const {
        Item item;
        pointToItem(q, 0, item);
        ConstIterator it =
            std::lower_bound(std::begin(items_), std::end(items_), item,
                             [&](const Item& item1, const Item& item2) -> bool {
                                 return mortonCmpInt(item1.index, item2.index);
                             });
        return it;
    }

    int computeLevel(ConstIterator first, ConstIterator last) {
        // For iterator syntax interval [first, last[ does not include last
        // iterator, but some operations requires access to last item. Hence, we
        // --last
        if (last != first) {
            --last;
        }
        return mortonDistance(first->index, last->index);
    }

    /**
     * Returns the iterator pointing to the splitting item of the given
     * input interval [first, last[.
     * In the implicit octree storage, there is an octree cell that contains
     * all the items of such interval.
     * Such cell is split into two octree cells [first, ret[ and [ret, last[
     * along the dimension where the original cell is larger.
     *
     * @param first the iterator pointing to the first item of the interval
     * @param last the iterator pointing after the last item of the interval
     * @param low the lower corner cell of the octree cell containing [first,
     * last[
     * @param mid the splitting cell (iterator points to the item NOT LESS than
     * this one)
     * @param upp the upper corner cell of the octree cell containing [first,
     * last[
     * @return
     */
    ConstIterator findSplit(ConstIterator first,
                            ConstIterator last,
                            PointInt& low,
                            PointInt& mid,
                            PointInt& upp) const {
        Item itemMid;
        if (last != first) {
            --last;
        }
        mortonSplit(first->index, last->index, low, mid, upp);
        itemMid.index = mid;
        auto it = std::lower_bound(items_.begin(), items_.end(), itemMid,
                                   [&](const Item& i1, const Item& i2) -> bool {
                                       return mortonCmpInt(i1.index, i2.index);
                                   });
        //            if (it != items_.end()) {
        //                ARS_VARIABLE(it->index.transpose());
        //            } else {
        //                ARS_PRINT("it == items_.end()");
        //            }
        return it;
    }

    /**
     * Returns the iterator pointing to the splitting item of the given
     * input interval [first, last[.
     * In the implicit octree storage, there is an octree cell that contains
     * all the items of such interval.
     * Such cell is split into two octree cells [first, ret[ and [ret, last[
     * along the dimension where the original cell is larger.
     *
     * @param first the iterator pointing to the first item of the interval
     * @param last the iterator pointing after the last item of the interval
     * @return
     */
    ConstIterator findSplit(ConstIterator first, ConstIterator last) const {
        PointInt low, mid, upp;
        return findSplit(first, last, low, mid, upp);
    }

    ConstIterator findSplit(const PointScalar& q0,
                            const PointScalar& q1) const {
        using MIT = MortonIntegerTraits<Integer>;
        using UnsignedType = typename MIT::UnsignedType;
        Item item0, item1, item01;
        PointInt low, upp;

        pointToItem(q0, 0, item0);
        pointToItem(q1, 0, item1);
        mortonSplit(item0.index, item1.index, low, item01, upp);

        return std::upper_bound(items_.begin(), items_.end(), item01,
                                [&](const Item& i1, const Item& i2) -> bool {
                                    return mortonCmpInt(i1.index, i2.index);
                                });
    }

   private:
    VectorItem items_;
    Scalar res_;

    void pointToItem(const PointScalar& p, size_t pos, Item& item) const {
        item.value = p;
        for (size_t d = 0; d < Dim; ++d) {
            item.index(d) = (Integer)round(p(d) / res_);
        }
        item.pos = pos;
    }
};

}  // namespace ars

#endif /* MORTONOCTREE_H */
