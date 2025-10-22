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
#include <map>

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
    using DiffT = std::map<int, double>;
    using FullT = std::vector<double>;
    using FullTPtr = std::shared_ptr<FullT>;

    ArsGraphInterval();

    ArsGraphInterval(ArsGraphPtr& graph);

    virtual ~ArsGraphInterval();

    size_t getNodeNum() const;

    size_t getEdgeNum() const;

    void setGraph(ArsGraphPtr& graph);

    virtual int getDiffSize() = 0;

    virtual double nodeLower(size_t i) = 0;

    virtual void setNodeLower(size_t idx, double val) = 0;

    virtual double nodeUpper(size_t i) = 0;

    virtual void setNodeUpper(size_t idx, double val) = 0;

    virtual std::vector<double> getNodeLowers() const = 0;

    virtual void setNodeLowers(std::vector<double>& vals) = 0;

    virtual void setNodeLowers(DiffT& vals) = 0;

    virtual std::vector<double> getNodeUppers() const = 0;

    virtual void setNodeUppers(std::vector<double>& vals) = 0;

    virtual void setNodeUppers(DiffT& vals) = 0;

    virtual double edgeLower(size_t i) = 0;

    virtual void setEdgeLower(size_t idx, double val) = 0;

    virtual double edgeUpper(size_t i) = 0;

    virtual void setEdgeUpper(size_t idx, double val) = 0;

    virtual std::vector<double> getEdgeLowers() const = 0;

    virtual void setEdgeLowers(std::vector<double>& vals) = 0;

    virtual void setEdgeLowers(DiffT& vals) = 0;

    virtual std::vector<double> getEdgeUppers() const = 0;

    virtual void setEdgeUppers(std::vector<double>& vals) = 0;

    virtual void setEdgeUppers(DiffT& vals) = 0;

    virtual double getLowerBound() = 0;

    virtual double getUpperBound() = 0;

    virtual void split(size_t idx, Ptr intervLower, Ptr intervUpper) = 0;

    virtual size_t size() const = 0;

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

    void init(ArsGraphPtr& graph);

    virtual int getDiffSize();

    virtual double nodeLower(size_t i);

    virtual void setNodeLower(size_t idx, double val);

    virtual double nodeUpper(size_t i);

    virtual void setNodeUpper(size_t idx, double val);

    virtual std::vector<double> getNodeLowers() const;

    virtual void setNodeLowers(std::vector<double>& vals);

    virtual void setNodeLowers(DiffT& vals);

    virtual std::vector<double> getNodeUppers() const;

    virtual void setNodeUppers(std::vector<double>& vals);

    virtual void setNodeUppers(DiffT& vals);

    virtual double edgeLower(size_t i);

    virtual void setEdgeLower(size_t idx, double val);

    virtual double edgeUpper(size_t i);

    virtual void setEdgeUpper(size_t idx, double val);

    virtual std::vector<double> getEdgeLowers() const;

    virtual void setEdgeLowers(std::vector<double>& vals);

    virtual void setEdgeLowers(DiffT& vals);

    virtual std::vector<double> getEdgeUppers() const;

    virtual void setEdgeUppers(std::vector<double>& vals);

    virtual void setEdgeUppers(DiffT& vals);

    virtual double getLowerBound();

    virtual double getUpperBound();

    virtual void getEdgeBounds(double& lower, double& upper);

    virtual void split(size_t idx,
                       Base::Ptr intervLower,
                       Base::Ptr intervUpper);
    
    virtual size_t size() const;    

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

class ArsGraphIntervalDiff : public ArsGraphInterval {
   public:
    using Base = ArsGraphInterval;
    using Self = ArsGraphIntervalDiff;
    using Ptr = std::shared_ptr<Self>;
    using ArsGraphPtr = ArsGraph::Ptr;

    ArsGraphIntervalDiff();

    ArsGraphIntervalDiff(ArsGraphPtr& graph);

    virtual ~ArsGraphIntervalDiff();

    virtual int getDiffSize();

    virtual double nodeLower(size_t i);

    virtual void setNodeLower(size_t idx, double val);

    virtual double nodeUpper(size_t i);

    virtual void setNodeUpper(size_t idx, double val);

    virtual std::vector<double> getNodeLowers() const;

    virtual void setNodeLowers(std::vector<double>& vals);

    virtual void setNodeLowers(DiffT& vals);

    virtual std::vector<double> getNodeUppers() const;

    virtual void setNodeUppers(std::vector<double>& vals);

    virtual void setNodeUppers(DiffT& vals);

    virtual double edgeLower(size_t i);

    virtual void setEdgeLower(size_t idx, double val);

    virtual double edgeUpper(size_t i);

    virtual void setEdgeUpper(size_t idx, double val);

    virtual std::vector<double> getEdgeLowers() const;

    virtual void setEdgeLowers(std::vector<double>& vals);

    virtual void setEdgeLowers(DiffT& vals);

    virtual std::vector<double> getEdgeUppers() const;

    virtual void setEdgeUppers(std::vector<double>& vals);

    virtual void setEdgeUppers(DiffT& vals);

    virtual double getLowerBound();

    virtual double getUpperBound();

    virtual void getEdgeBounds(double& lower, double& upper);

    virtual void split(size_t idx,
                       Base::Ptr intervLower,
                       Base::Ptr intervUpper);
    
    virtual size_t size() const;    

   private:
    FullTPtr nodeLowersParent_;
    DiffT nodeLowersDiff_;

    FullTPtr nodeUppersParent_;
    DiffT nodeUppersDiff_;

    FullTPtr edgeLowersParent_;
    DiffT edgeLowersDiff_;

    FullTPtr edgeUppersParent_;
    DiffT edgeUppersDiff_;

    double lower_;
    double upper_;

    void computeEdgeBounds();
};

}  // namespace ars

#endif