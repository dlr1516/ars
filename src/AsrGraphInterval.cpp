#include <ars/ArsGraphInterval.h>

#include < cmath>

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
    return graph_->getNodeNum();
}

// ------------------------------------------------------------------
// ARS GRAPH INTERVALL FULL
// ------------------------------------------------------------------

ArsGraphIntervalFull::ArsGraphIntervalFull() : ArsGraphInterval() {}

ArsGraphIntervalFull::ArsGraphIntervalFull(ArsGraphPtr& graph)
    : ArsGraphInterval(graph_) {}

ArsGraphIntervalFull::~ArsGraphIntervalFull() {}

void ArsGraphIntervalFull::init(ArsGraphPtr& graph) {
    graph_ = graph;

    size_t nodeNum = this->getNodeNum();
    size_t edgeNum = this->getEdgeNum();

    // Initializes node lower and upper values
    nodeLowers_.resize(nodeNum);
    nodeUppers_.resize(nodeNum);
    std::fill(std::begin(nodeLowers_), std::end(nodeLowers_), 0.0);
    std::fill(std::begin(nodeUppers_), std::end(nodeUppers_), M_PI);

    // Computes lower and upper bounds of each edge
    edgeLowers_.resize(edgeNum);
    edgeUppers_.resize(edgeNum);
    for (int j = 0; j < graph_->edges().size(); ++j) {
        ArsGraph::Edge& edge = graph_->edges()[j];
        int isrc = edge.isrc;
        int idst = edge.idst;
        double thetaMin = nodeLowers_[idst] - nodeUppers_[isrc];
        double thetaMax = nodeUppers_[idst] - nodeLowers_[isrc];
        ars::findLUFourier(edge.coeffs, thetaMin, thetaMax, , edgeLowers_[j],
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
                                 ArsGraphInterval::Ptr intervLower,
                                 ArsGraphInterval::Ptr intervUpper) {
    intervLower = std::make_shared<ArsGraphIntervalFull>();
    intervUpper = std::make_shared<ArsGraphIntervalFull>();

    intervLower->nodeLowers_ = nodeLowers_;
    intervLower->nodeUppers_ = nodeUppers_;

    intervLower->nodeLowers_ = nodeLowers_;
    intervLower->nodeUppers_ = nodeUppers_;
}

};  // namespace ars