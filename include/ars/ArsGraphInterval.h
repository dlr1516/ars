/**
 * ARS - Angular Radon Spectrum
 * Copyright (C) 2025 Dario Lodi Rizzini - Ernesto Fontana.
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
#ifndef ARSGRAPHINTERVAL_H
#define ARSGRAPHINTERVAL_H

#include <memory>
#include <vector>

#include <ars/ArsGraph.h>
#include <ars/definitions.h>
#include <ars/functions.h>

namespace ars {

// --------------------------------------------------------
// ARS GRAPH INTERVAL
// --------------------------------------------------------

class ArsGraphInterval {
   public:
    using Self = ArsGraphInterval;
    using Ptr = std::shared_ptr<Self>;
    using ArsGraphPtr = typename ArsGraph::Ptr;

    ArsGraphInterval();

    ArsGraphInterval(ArsGraphPtr& graph);

    virtual ~ArsGraphInterval();

    size_t getNodeNum() const;

    size_t getEdgeNum() const;

    virtual double stateLower(size_t i) = 0;

    virtual double stateUpper(size_t i) = 0;

    void split(size_t i, Ptr intervLow, Ptr intervUpp);

   private:
    ArsGraphPtr graph_;
};

// --------------------------------------------------------
// ARS GRAPH INTERVAL FULL
// --------------------------------------------------------

class ArsGraphIntervalFull : public ArsGraphInterval {
   public:
    using Base = ArsGraphInterval;
    using Self = ArsGraphIntervalFull;
    using Ptr = std::shared_ptr<Self>;
    using ArsGraphPtr = typename ArsGraph::Ptr;

    ArsGraphIntervalFull();

    ArsGraphIntervalFull(ArsGraphPtr& graph);

    ~ArsGraphIntervalFull();

    void init(ArsGraphPtr& graph);

    double stateLower(size_t i);

    double stateUpper(size_t i);

    void getEdgeBounds(double& lower, double& upper);

    void split(int idx,
                       ArsGraphIntervalFull::Ptr intervLower,
                       ArsGraphIntervalFull::Ptr intervUpper);

   private:
    ArsGraphPtr graph_;
    std::vector<double> nodeLowers_;
    std::vector<double> nodeUppers_;
    std::vector<double> edgeLowers_;
    std::vector<double> edgeUppers_;
};

// --------------------------------------------------------
// ARS GRAPH INTERVAL DIFF
// --------------------------------------------------------

}  // namespace ars

#endif