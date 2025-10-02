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
    using ArsGraphPtr = ArsGraph::Ptr;

    ArsGraphInterval();

    ArsGraphInterval(ArsGraphPtr& graph);

    virtual ~ArsGraphInterval();

    size_t getNodeNum() const;

    size_t getEdgeNum() const;

    virtual double stateLower(size_t i) = 0;

    virtual double stateUpper(size_t i) = 0;

    virtual std::vector<double> getNodesLower() const = 0;

    virtual std::vector<double> getNodesUpper() const = 0;

    virtual double getLower() const = 0;

    virtual double getUpper() const = 0;

    void split(size_t i, Ptr intervLow, Ptr intervUpp);

   protected:
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
    using ArsGraphPtr = ArsGraph::Ptr;

    ArsGraphIntervalFull();

    ArsGraphIntervalFull(ArsGraphPtr& graph);

    virtual ~ArsGraphIntervalFull();

    virtual void init(ArsGraphPtr& graph);

    virtual double stateLower(size_t i);

    virtual double stateUpper(size_t i);

    virtual std::vector<double> getNodesLower() const;

    virtual std::vector<double> getNodesUpper() const;

    virtual double edgeLower(size_t i);

    virtual double edgeUpper(size_t i);

    virtual double getLower() const;

    virtual double getUpper() const;

    virtual void getEdgeBounds(double& lower, double& upper);

    virtual void split(int idx,
                       ArsGraphIntervalFull::Ptr& intervLower,
                       ArsGraphIntervalFull::Ptr& intervUpper);

   private:
    std::vector<double> nodeLowers_;
    std::vector<double> nodeUppers_;
    std::vector<double> edgeLowers_;
    std::vector<double> edgeUppers_;
    double lower_;
    double upper_;
};

// --------------------------------------------------------
// ARS GRAPH INTERVAL DIFF
// --------------------------------------------------------

}  // namespace ars

#endif