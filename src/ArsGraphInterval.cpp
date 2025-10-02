#include <ars/ArsGraphInterval.h>

#include <cmath>

namespace ars {

// ------------------------------------------------------------------
// ARS GRAPH INTERVALL BASE
// ------------------------------------------------------------------

ArsGraphInterval::ArsGraphInterval() : graph_(nullptr) {}

ArsGraphInterval::ArsGraphInterval(ArsGraphPtr& graph) : graph_(graph) {}

ArsGraphInterval::~ArsGraphInterval() {}

size_t ArsGraphInterval::getNodeNum() const {
    if (graph_.get() == nullptr) {
        ARS_ERROR("interval not initialized");
        return 0;
    }
    return graph_->getNodeNum();
}

size_t ArsGraphInterval::getEdgeNum() const {
    if (graph_.get() == nullptr) {
        ARS_ERROR("interval not initialized");
        return 0;
    }
    return graph_->getEdgeNum();
}

// ------------------------------------------------------------------
// ARS GRAPH INTERVALL FULL
// ------------------------------------------------------------------

ArsGraphIntervalFull::ArsGraphIntervalFull() : ArsGraphInterval() {}

ArsGraphIntervalFull::ArsGraphIntervalFull(ArsGraphPtr& graph)
    : ArsGraphInterval(graph) {}

ArsGraphIntervalFull::~ArsGraphIntervalFull() {}

void ArsGraphIntervalFull::init(ArsGraphPtr& graph) {
    graph_ = graph;

    size_t nodeNum = getNodeNum();
    size_t edgeNum = getEdgeNum();

    // Initializes node lower and upper values
    nodeLowers_.resize(nodeNum);
    nodeUppers_.resize(nodeNum);
    std::fill(std::begin(nodeLowers_), std::end(nodeLowers_), 0.0);
    std::fill(std::begin(nodeUppers_), std::end(nodeUppers_), M_PI);
    nodeLowers_[0] = 0.0;
    nodeUppers_[0] = 0.0;

    // Computes lower and upper bounds of each edge
    edgeLowers_.resize(edgeNum);
    edgeUppers_.resize(edgeNum);
    for (int j = 0; j < graph_->edges().size(); ++j) {
        const ArsGraph::Edge& edge = graph_->edges()[j];
        int isrc = edge.isrc;
        int idst = edge.idst;
        double thetaMin = nodeLowers_[idst] - nodeUppers_[isrc];
        double thetaMax = nodeUppers_[idst] - nodeLowers_[isrc];
        ars::findLUFourier(edge.coeffs, thetaMin, thetaMax, edgeLowers_[j],
                           edgeUppers_[j]);
    }
}

double ArsGraphIntervalFull::stateLower(size_t i) {
    // TODO: check overflow of i
    return nodeLowers_[i];
}

double ArsGraphIntervalFull::stateUpper(size_t i) {
    // TODO: check overflow of i
    return nodeUppers_[i];
}

std::vector<double> ArsGraphIntervalFull::getNodesLower() const {
    return nodeLowers_;
}

std::vector<double> ArsGraphIntervalFull::getNodesUpper() const {
    return nodeUppers_;
}

double ArsGraphIntervalFull::edgeLower(size_t i) {
    // TODO: check overflow of i
    return edgeLowers_[i];
}

double ArsGraphIntervalFull::edgeUpper(size_t i) {
    // TODO: check overflow of i
    return edgeUppers_[i];
}

double ArsGraphIntervalFull::getLower() const {
    return lower_;
}

double ArsGraphIntervalFull::getUpper() const {
    return lower_;
}

void ArsGraphIntervalFull::getEdgeBounds(double& lower, double& upper) {
    int edgeNum = edgeLowers_.size();
    lower = 0.0;
    upper = 0.0;
    for (int j = 0; j < edgeNum; ++j) {
        lower += edgeLowers_[j];
        upper += edgeUppers_[j];
    }
}

void ArsGraphIntervalFull::split(int idx,
                                 ArsGraphIntervalFull::Ptr& intervLower,
                                 ArsGraphIntervalFull::Ptr& intervUpper) {
    ARS_ASSERT_VAR1(idx != 0, idx);

    intervLower =
        std::make_shared<ArsGraphIntervalFull>(ArsGraphIntervalFull(graph_));
    intervUpper =
        std::make_shared<ArsGraphIntervalFull>(ArsGraphIntervalFull(graph_));

    double thetaMid = 0.5 * (nodeLowers_[idx] + nodeUppers_[idx]);
    intervLower->nodeLowers_ = nodeLowers_;
    intervLower->nodeUppers_ = nodeUppers_;
    intervLower->nodeUppers_[idx] = thetaMid;
    intervUpper->nodeLowers_ = nodeLowers_;
    intervUpper->nodeUppers_ = nodeUppers_;
    intervUpper->nodeLowers_[idx] = thetaMid;

    intervLower->edgeLowers_ = edgeLowers_;
    intervLower->edgeUppers_ = edgeUppers_;
    intervUpper->edgeLowers_ = edgeLowers_;
    intervUpper->edgeUppers_ = edgeUppers_;

    for (int eidx : graph_->nodes()[idx].incidents) {
        const ArsGraph::Edge& edge = graph_->edges()[eidx];
        int isrc = edge.isrc;
        int idst = edge.idst;

        double thetaMin =
            intervLower->nodeLowers_[idst] - intervLower->nodeUppers_[isrc];
        double thetaMax =
            intervLower->nodeUppers_[idst] - intervLower->nodeLowers_[isrc];
        ars::findLUFourier(edge.coeffs, thetaMin, thetaMax,
                           intervLower->edgeLowers_[eidx],
                           intervLower->edgeUppers_[eidx]);

        thetaMin =
            intervUpper->nodeLowers_[idst] - intervUpper->nodeUppers_[isrc];
        thetaMax =
            intervUpper->nodeUppers_[idst] - intervUpper->nodeLowers_[isrc];
        ars::findLUFourier(edge.coeffs, thetaMin, thetaMax,
                           intervUpper->edgeLowers_[eidx],
                           intervUpper->edgeUppers_[eidx]);
    }

    intervLower->getEdgeBounds(intervLower->lower_, intervLower->upper_);
    intervUpper->getEdgeBounds(intervUpper->lower_, intervUpper->upper_);
}

};  // namespace ars