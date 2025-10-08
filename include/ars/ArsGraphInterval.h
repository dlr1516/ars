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

    void setGraph(ArsGraphPtr& graph);

    virtual double nodeLower(size_t i) = 0;

    virtual void setNodeLower(size_t idx, double val) = 0;

    virtual double nodeUpper(size_t i) = 0;

    virtual void setNodeUpper(size_t idx, double val) = 0;

    virtual std::vector<double> getNodesLower() const = 0;

    virtual void setNodesLower(std::vector<double> vals) = 0;

    virtual std::vector<double> getNodesUpper() const = 0;

    virtual void setNodesUpper(std::vector<double> vals) = 0;

    virtual double edgeLower(size_t i) = 0;

    virtual void setEdgeLower(size_t idx, double val) = 0;

    virtual double edgeUpper(size_t i) = 0;

    virtual void setEdgeUpper(size_t idx, double val) = 0;

    virtual std::vector<double> getEdgesLower() const = 0;

    virtual void setEdgesLower(std::vector<double> vals) = 0;

    virtual std::vector<double> getEdgesUpper() const = 0;

    virtual void setEdgesUpper(std::vector<double> vals) = 0;

    virtual double getLowerBound() = 0;

    virtual double getUpperBound() = 0;

    virtual void split(size_t idx, Ptr intervLower, Ptr intervUpper) = 0;

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

    virtual double nodeLower(size_t i);

    virtual void setNodeLower(size_t idx, double val);

    virtual double nodeUpper(size_t i);

    virtual void setNodeUpper(size_t idx, double val);

    virtual std::vector<double> getNodesLower() const;

    virtual void setNodesLower(std::vector<double> vals);

    virtual std::vector<double> getNodesUpper() const;

    virtual void setNodesUpper(std::vector<double> vals);

    virtual double edgeLower(size_t i);

    virtual void setEdgeLower(size_t idx, double val);

    virtual double edgeUpper(size_t i);

    virtual void setEdgeUpper(size_t idx, double val);

    virtual std::vector<double> getEdgesLower() const;

    virtual void setEdgesLower(std::vector<double> vals);

    virtual std::vector<double> getEdgesUpper() const;

    virtual void setEdgesUpper(std::vector<double> vals);

    virtual double getLowerBound();

    virtual double getUpperBound();

    virtual void getEdgeBounds(double& lower, double& upper);

    virtual void split(size_t idx,
                       ArsGraphInterval::Ptr intervLower,
                       ArsGraphInterval::Ptr intervUpper);

   private:
    std::vector<double> nodeLowers_;
    std::vector<double> nodeUppers_;
    std::vector<double> edgeLowers_;
    std::vector<double> edgeUppers_;
    double lower_;
    double upper_;

    void computeEdgeBounds();
};

// --------------------------------------------------------
// ARS GRAPH INTERVAL DIFF
// --------------------------------------------------------

}  // namespace ars

#endif