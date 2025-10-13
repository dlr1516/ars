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

void ArsGraphInterval::setGraph(ArsGraphPtr& graph){
    graph_ = graph;
}

// ------------------------------------------------------------------
// ARS GRAPH INTERVALL FULL
// ------------------------------------------------------------------

ArsGraphIntervalFull::ArsGraphIntervalFull() 
    : ArsGraphInterval(), upper_(NAN), lower_(NAN) {}

ArsGraphIntervalFull::ArsGraphIntervalFull(ArsGraphPtr& graph)
    : ArsGraphInterval(graph), upper_(NAN), lower_(NAN) {}

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
    computeEdgeBounds();
}

double ArsGraphIntervalFull::nodeLower(size_t i) {
    // TODO: check overflow of i
    return nodeLowers_[i];
}

void ArsGraphIntervalFull::setNodeLower(size_t idx, double val) {
    nodeLowers_[idx] = val;
}

double ArsGraphIntervalFull::nodeUpper(size_t i)
{
    // TODO: check overflow of i
    return nodeUppers_[i];
}

void ArsGraphIntervalFull::setNodeUpper(size_t idx, double val) {
    nodeUppers_[idx] = val;
}

std::vector<double> ArsGraphIntervalFull::getNodeLowers() const
{
    return nodeLowers_;
}

void ArsGraphIntervalFull::setNodeLowers(std::vector<double>& vals) {
    nodeLowers_ = vals;
}

std::vector<double> ArsGraphIntervalFull::getNodeUppers() const {
    return nodeUppers_;
}

void ArsGraphIntervalFull::setNodeUppers(std::vector<double>& vals) {
    nodeUppers_ = vals;
}

double ArsGraphIntervalFull::edgeLower(size_t i) {
    // TODO: check overflow of i
    return edgeLowers_[i];
}

void ArsGraphIntervalFull::setEdgeLower(size_t idx, double val) {
    edgeLowers_[idx] = val;
}

double ArsGraphIntervalFull::edgeUpper(size_t i) {
    // TODO: check overflow of i
    return edgeUppers_[i];
}

void ArsGraphIntervalFull::setEdgeUpper(size_t idx, double val) {
    edgeUppers_[idx] = val;
}

std::vector<double> ArsGraphIntervalFull::getEdgeLowers() const {
    return edgeLowers_;
}

void ArsGraphIntervalFull::setEdgeLowers(std::vector<double>& vals) {
    edgeLowers_ = vals;
}

std::vector<double> ArsGraphIntervalFull::getEdgeUppers() const {
    return edgeUppers_;
}

void ArsGraphIntervalFull::setEdgeUppers(std::vector<double>& vals) {
    edgeUppers_ = vals;
}

double ArsGraphIntervalFull::getLowerBound() {
    if(std::isnan(lower_)){
        computeEdgeBounds();
    }
    return lower_;
}

double ArsGraphIntervalFull::getUpperBound() {
    if(std::isnan(upper_)){
        computeEdgeBounds();
    }
    return upper_;
}

void ArsGraphIntervalFull::getEdgeBounds(double& lower, double& upper) {
    if(std::isnan(lower_) || std::isnan(upper_)){
        computeEdgeBounds();
    }
    lower = lower_;
    upper = upper_;
}

void ArsGraphIntervalFull::computeEdgeBounds() {
    int edgeNum = edgeLowers_.size();
    lower_ = 0.0;
    upper_ = 0.0;
    for (int i = 0; i < edgeNum; ++i) {
        lower_ += edgeLowers_[i];
        upper_ += edgeUppers_[i];
    }
}

void ArsGraphIntervalFull::split(size_t idx,
                                 ArsGraphInterval::Ptr intervLower,
                                 ArsGraphInterval::Ptr intervUpper) {
    ARS_ASSERT_VAR1(idx != 0, idx);

    intervLower->setGraph(graph_);
    intervUpper->setGraph(graph_);

    double thetaMid = 0.5 * (nodeLowers_[idx] + nodeUppers_[idx]);
    intervLower->setNodeLowers(nodeLowers_);
    intervLower->setNodeUppers(nodeUppers_);
    intervLower->setNodeUpper(idx, thetaMid);

    intervUpper->setNodeLowers(nodeLowers_);
    intervUpper->setNodeUppers(nodeUppers_);
    intervUpper->setNodeLower(idx, thetaMid);

    intervLower->setEdgeLowers(edgeLowers_);
    intervLower->setEdgeUppers(edgeUppers_);
    intervUpper->setEdgeLowers(edgeLowers_);
    intervUpper->setEdgeUppers(edgeUppers_);

    for (int eidx : graph_->nodes()[idx].incidents) {
        const ArsGraph::Edge& edge = graph_->edges()[eidx];
        int isrc = edge.isrc;
        int idst = edge.idst;
        double edgeLower, edgeUpper;

        double thetaMin =
            intervLower->nodeLower(idst) - intervLower->nodeUpper(isrc);
        double thetaMax =
            intervLower->nodeUpper(idst) - intervLower->nodeLower(isrc);
        ars::findLUFourier(edge.coeffs, thetaMin, thetaMax,
                           edgeLower, edgeUpper);
        intervLower->setEdgeLower(eidx, edgeLower);
        intervLower->setEdgeUpper(eidx, edgeUpper);

        thetaMin =
            intervUpper->nodeLower(idst) - intervUpper->nodeUpper(isrc);
        thetaMax =
            intervUpper->nodeUpper(idst) - intervUpper->nodeLower(isrc);
        ars::findLUFourier(edge.coeffs, thetaMin, thetaMax,
                           edgeLower, edgeUpper);
        intervUpper->setEdgeLower(eidx, edgeLower);
        intervUpper->setEdgeUpper(eidx, edgeUpper);
    }
}

size_t ArsGraphIntervalFull::size() const{
    size_t size = 0;
    size += sizeof(graph_);
    size += nodeLowers_.size()*sizeof(nodeLowers_.front());
    size += nodeUppers_.size()*sizeof(nodeUppers_.front());
    size += edgeLowers_.size()*sizeof(edgeLowers_.front());
    size += edgeUppers_.size()*sizeof(edgeUppers_.front());
    size += sizeof(lower_);
    size += sizeof(upper_);
    return size;
}
}; // namespace ars