/**
 * ARS - Angular Radon Spectrum
 * Copyright (C) 2017 Dario Lodi Rizzini.
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
#ifndef CONSENSUS_TRANSLATION_ESTIMATOR_H
#define CONSENSUS_TRANSLATION_ESTIMATOR_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <ars/definitions.h>
#include <rofl/common/grid.h>
#include <rofl/common/peak_finder_d.h>

#include <boost/range/adaptor/indexed.hpp>
using namespace boost::adaptors;

namespace ars
{

    /**
     * Class ConsensusTranslationEstimator provides a simple solution for
     * computing the translation between two point sets, which we assume to
     * have the same orientation.
     * The principle is simple: given the source points pointsSrc[i],
     * the destination points pointsDst[i], with perfect correspondence between them,
     * they are related by translation t
     *
     *  pointsDst[i] = pointsSrc[i] + t  -> t = pointsDst[i] - pointsSrc[i]
     *
     * Since we do not have such clear correspondence (some outlier may not have
     * association), what we compute all the vectors
     *
     *   t[i][j] = pointsDst[j] - pointsSrc[i]
     *
     * When the pair (i, j) are right correspondences, t[i][j] is (close to)
     * the translation value. We assume that this is the majority of cases.
     * Thus, ConsensusTranslationEstimator computes the consensus using a
     * classic voting grid.
     * Translation vector candidates correspond to maxima in grid.
     */
    template <size_t Dim, typename Scalar = double>
    class ConsensusTranslationEstimator
    {
    public:
        using Index = int;
        using Counter = size_t;
        using Grid = rofl::Grid<Dim, Counter, Index, rofl::detail::RasterIndexer<2, Index>, std::vector, std::allocator>;
        using Indices = typename Grid::Indices;
        using PeakFinder = rofl::PeakFinderD<Dim, Counter, Index, std::greater<Index>>;

#if __cplusplus < 201703L
        using Point = Eigen::Matrix<Scalar, Dim, 1>;
        using VectorPoint = std::vector<Point, Eigen::aligned_allocator<Point>>;
#else
        using Point = Eigen::Matrix<Scalar, Dim, 1>;
        using VectorPoint = std::vector<Point>;
#endif

        /**
         * Default constructor.
         */
        ConsensusTranslationEstimator() : grid_(), gridIndices_(Dim), translMin_(Point::Zero()), translRes_(1.0), peakFinder_()
        {
        }

        ConsensusTranslationEstimator(const Point &translMin, const Scalar &translRes, const Indices &gridSize)
            : grid_(), gridIndices_(Dim), translMin_(translMin), translRes_(translRes), peakFinder_()
        {
            grid_.initBounds(gridSize);
            peakFinder_.setDomain(gridSize);
        }

        virtual ~ConsensusTranslationEstimator()
        {
        }

        void init(const Point &translMin, const Scalar &translRes, const Indices &gridSize)
        {
            grid_.initBounds(gridSize);
            translMin_ = translMin;
            translRes_ = translRes;
            peakFinder_.setDomain(gridSize);
        }

        void reset()
        {
            grid_.fill(0);

            if (gridIndices_.empty())
                return;
            for (int i = 0; i < Dim; ++i)
                if (!gridIndices_[i].empty())
                    gridIndices_[i].clear();
            // gridIndices_.clear();
        }

        void setNonMaximaWindowDim(const Indices &dim)
        {
            peakFinder_.setPeakWindow(dim);
        }

        void insert(const VectorPoint &pointsSrc, const VectorPoint &pointsDst, bool adaptive = false)
        {
            Point translMax, srcMin, srcMax, dstMin, dstMax;
            Indices indices, gridSize;

            if (adaptive)
            {
                srcMin.fill(std::numeric_limits<Scalar>::max());
                srcMax.fill(std::numeric_limits<Scalar>::lowest());
                for (auto &p : pointsSrc)
                {
                    for (int d = 0; d < Dim; ++d)
                    {
                        if (p(d) < srcMin(d))
                            srcMin(d) = p(d);
                        if (p(d) > srcMax(d))
                            srcMax(d) = p(d);
                    }
                }
                dstMin.fill(std::numeric_limits<Scalar>::max());
                dstMax.fill(std::numeric_limits<Scalar>::lowest());
                for (auto &p : pointsDst)
                {
                    for (int d = 0; d < Dim; ++d)
                    {
                        if (p(d) < dstMin(d))
                            dstMin(d) = p(d);
                        if (p(d) > dstMax(d))
                            dstMax(d) = p(d);
                    }
                }
                translMin_ = dstMin - srcMax;
                translMax = dstMax - srcMin;
                for (int d = 0; d < Dim; ++d)
                {
                    gridSize[d] = (Index)ceil((translMax(d) - translMin_(d)) / translRes_);
                }
                //                ARS_VAR5(translMin_.transpose(), translMax.transpose(), translRes_, gridSize[0], gridSize[1]);
                init(translMin_, translRes_, gridSize);
            }

            gridIndices_.resize(grid_.size()); // FIXME!

            int counterMax = 0;
            Indices indicesMax;
            for (int i = 0; i < pointsSrc.size(); ++i)
            {
                for (int j = 0; j < pointsDst.size(); ++j)
                {
                    Point ps = pointsSrc[i];
                    Point pd = pointsDst[i];
                    Point transl = pd - ps;
                    indices = getIndices(transl);
                    // ARS_VARIABLE4(transl.transpose(),indices[0],indices[1],grid_.inside(indices));
                    if (grid_.inside(indices))
                    {
                        grid_.value(indices)++;
                        if (grid_.value(indices) > counterMax)
                        {
                            counterMax = grid_.value(indices);
                            indicesMax = indices;
                        }

                        gridIndices_[grid_.getPos(indices)].push_back(std::pair<int, int>(i, j));
                    }
                }
            }
            ARS_VAR3(counterMax, indicesMax[0], indicesMax[1]);
            Counter thres = std::min(pointsSrc.size(), pointsDst.size()) / 2;
            peakFinder_.enableFilterPeakMin(true, thres);
        }

        void computeMaxima(std::vector<Indices> &indicesMax)
        {
            indicesMax.clear();
            auto histoMap = [&](const Indices &indices) -> Counter
            {
                // ARS_VARIABLE3(indices[0], indices[1], grid_.inside(indices));
                return grid_.value(indices);
            };
            peakFinder_.detect(histoMap, std::back_inserter(indicesMax));
            //            ARS_PRINT("Maxima:");
            //            for (auto &idx : indicesMax) {
            //                std::cout << "  indices [" << idx[0] << "," << idx[1] << "] value " << histoMap(idx)
            //                        << " grid2.value() " << grid_.value(idx) << std::endl;
            //            }
        }

        void computeMaxima(VectorPoint &translMax, std::vector<Indices> &indicesMax)
        {
            computeMaxima(indicesMax);
            translMax.clear();
            translMax.reserve(indicesMax.size());
            for (auto idx : indicesMax)
            {
                Point p = getTranslation(idx);
                //                ARS_VAR4(idx[0], idx[1], grid_.value(idx), p.transpose());
                translMax.push_back(p);
            }
        }

        Indices getIndices(const Point &p) const
        {
            Indices indices;
            for (int d = 0; d < Dim; ++d)
            {
                indices[d] = round((p(d) - translMin_(d)) / translRes_);
            }
            return indices;
        }

        Point getTranslation(const Indices &indices) const
        {
            Point transl;
            for (int d = 0; d < Dim; ++d)
            {
                transl(d) = translRes_ * indices[d] + translMin_(d);
            }
            return transl;
        }

        Counter getScore(const Point &p) const
        {
            Indices indices = getIndices(p);
            return getScore(indices);
        }

        Counter getScore(const Indices &indices) const
        {
            return grid_.value(indices);
        }

        const Grid &getGrid() const
        {
            return grid_;
        }

        const std::vector<std::vector<std::pair<Index, Index>>> &getGridIndices() const
        {
            return gridIndices_;
        }

    private:
        Grid grid_;
        std::vector<std::vector<std::pair<Index, Index>>> gridIndices_; // TODO: make gridIndices_ explicitly of type rofl::Grid<pair>
        Point translMin_;
        Scalar translRes_;
        PeakFinder peakFinder_;
    };

    using ConsensusTranslationEstimator2f = ConsensusTranslationEstimator<2, float>;
    using ConsensusTranslationEstimator2d = ConsensusTranslationEstimator<2, double>;
    using ConsensusTranslationEstimator3f = ConsensusTranslationEstimator<3, float>;
    using ConsensusTranslationEstimator3d = ConsensusTranslationEstimator<3, double>;

} // end of namespace

#endif /* VOTINGTRANSLATIONESTIMATOR_H */
