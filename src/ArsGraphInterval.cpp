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

int ArsGraphIntervalFull::getDiffSize() {
    return 0;
}

double ArsGraphIntervalFull::nodeLower(size_t i)
{
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

void ArsGraphIntervalFull::setNodeLowers(std::map<int, double> &vals) {
    for (auto& [key, val] : vals)
        nodeLowers_[key] = val;
}

std::vector<double> ArsGraphIntervalFull::getNodeUppers() const
{
    return nodeUppers_;
}

void ArsGraphIntervalFull::setNodeUppers(std::vector<double>& vals) {
    nodeUppers_ = vals;
}

void ArsGraphIntervalFull::setNodeUppers(std::map<int, double> &vals) {
    for (auto& [key, val] : vals)
        nodeUppers_[key] = val;
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

void ArsGraphIntervalFull::setEdgeLowers(std::map<int, double> &vals) {
    for (auto& [key, val] : vals)
        edgeLowers_[key] = val;
}

std::vector<double> ArsGraphIntervalFull::getEdgeUppers() const {
    return edgeUppers_;
}

void ArsGraphIntervalFull::setEdgeUppers(std::vector<double>& vals) {
    edgeUppers_ = vals;
}

void ArsGraphIntervalFull::setEdgeUppers(std::map<int, double> &vals) {
    for (auto& [key, val] : vals)
        edgeUppers_[key] = val;
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
    int edgeNum = getEdgeNum();
    lower_ = 0.0;
    upper_ = 0.0;
    for (int i = 0; i < edgeNum; ++i) {
        lower_ += edgeLowers_[i];
        upper_ += edgeUppers_[i];
    }
}

void ArsGraphIntervalFull::split(size_t idx,
                                 Base::Ptr intervLower,
                                 Base::Ptr intervUpper) {
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

// ------------------------------------------------------------------
// ARS GRAPH INTERVALL DIFF
// ------------------------------------------------------------------

ArsGraphIntervalDiff::ArsGraphIntervalDiff() 
    : ArsGraphInterval(), upper_(NAN), lower_(NAN) {}

ArsGraphIntervalDiff::ArsGraphIntervalDiff(ArsGraphPtr& graph)
    : ArsGraphInterval(graph), upper_(NAN), lower_(NAN) {}

ArsGraphIntervalDiff::~ArsGraphIntervalDiff() {}

int ArsGraphIntervalDiff::getDiffSize() {
    return edgeLowersDiff_.size();
}

double ArsGraphIntervalDiff::nodeLower(size_t i) {
    if (nodeLowersDiff_.contains(i))
        return nodeLowersDiff_[i];
    else
        return nodeLowersParent_->at(i);
}

void ArsGraphIntervalDiff::setNodeLower(size_t idx, double val) {
    nodeLowersDiff_[idx] = val;
}

double ArsGraphIntervalDiff::nodeUpper(size_t i) {
    if (nodeUppersDiff_.contains(i))
        return nodeUppersDiff_[i];
    else
        return nodeUppersParent_->at(i);
}

void ArsGraphIntervalDiff::setNodeUpper(size_t idx, double val) {
    nodeUppersDiff_[idx] = val;
}

std::vector<double> ArsGraphIntervalDiff::getNodeLowers() const {
    std::vector<double> nodes = *nodeLowersParent_;
    for (auto& [key, val] : nodeLowersDiff_)
        nodes[key] = val;
    return nodes;
}

void ArsGraphIntervalDiff::setNodeLowers(std::vector<double> &vals) {
    nodeLowersParent_ = std::make_shared<FullT>(vals);
    nodeLowersDiff_.clear();
}

void ArsGraphIntervalDiff::setNodeLowers(std::map<int, double> &vals) {
    for (auto& [key, val] : vals)
        nodeLowersDiff_[key] = val;
}

std::vector<double> ArsGraphIntervalDiff::getNodeUppers() const {
    std::vector<double> nodes = *nodeUppersParent_;
    for (auto& [key, val] : nodeUppersDiff_)
        nodes[key] = val;
    return nodes;
}

void ArsGraphIntervalDiff::setNodeUppers(std::vector<double> &vals) {
    nodeUppersParent_ = std::make_shared<FullT>(vals);
    nodeUppersDiff_.clear();
}

void ArsGraphIntervalDiff::setNodeUppers(std::map<int, double> &vals) {
    for (auto& [key, val] : vals)
        nodeUppersDiff_[key] = val;
}

double ArsGraphIntervalDiff::edgeLower(size_t i)
{
    if (edgeLowersDiff_.contains(i))
        return edgeLowersDiff_[i];
    else
        return edgeLowersParent_->at(i);
}

void ArsGraphIntervalDiff::setEdgeLower(size_t idx, double val) {
    edgeLowersDiff_[idx] = val;
}

double ArsGraphIntervalDiff::edgeUpper(size_t i) {
    if (edgeUppersDiff_.contains(i))
        return edgeUppersDiff_[i];
    else
        return edgeUppersParent_->at(i);
}

void ArsGraphIntervalDiff::setEdgeUpper(size_t idx, double val) {
    edgeUppersDiff_[idx] = val;
}

std::vector<double> ArsGraphIntervalDiff::getEdgeLowers() const {
    std::vector<double> edges = *edgeLowersParent_;
    for (auto& [key, val] : edgeLowersDiff_)
        edges[key] = val;
    return edges;
}

void ArsGraphIntervalDiff::setEdgeLowers(std::vector<double> &vals) {
    edgeLowersParent_ = std::make_shared<FullT>(vals);
    edgeLowersDiff_.clear();
}

void ArsGraphIntervalDiff::setEdgeLowers(std::map<int, double> &vals) {
    for (auto& [key, val] : vals)
        edgeLowersDiff_[key] = val;
}

std::vector<double> ArsGraphIntervalDiff::getEdgeUppers() const {
    std::vector<double> edges = *edgeUppersParent_;
    for (auto& [key, val] : edgeUppersDiff_)
        edges[key] = val;
    return edges;
}

void ArsGraphIntervalDiff::setEdgeUppers(std::vector<double> &vals) {
    edgeUppersParent_ = std::make_shared<FullT>(vals);
    edgeUppersDiff_.clear();
}

void ArsGraphIntervalDiff::setEdgeUppers(std::map<int, double> &vals) {
    for (auto& [key, val] : vals)
        edgeUppersDiff_[key] = val;
}

double ArsGraphIntervalDiff::getLowerBound() {
    if(std::isnan(lower_)){
        computeEdgeBounds();
    }
    return lower_;
}

double ArsGraphIntervalDiff::getUpperBound() {
    if(std::isnan(upper_)){
        computeEdgeBounds();
    }
    return upper_;
}

void ArsGraphIntervalDiff::getEdgeBounds(double& lower, double& upper) {
    if(std::isnan(lower_) || std::isnan(upper_)){
        computeEdgeBounds();
    }
    lower = lower_;
    upper = upper_;
}

void ArsGraphIntervalDiff::computeEdgeBounds() {
    int edgeNum = getEdgeNum();
    lower_ = 0.0;
    upper_ = 0.0;
    for (int i = 0; i < edgeNum; ++i) {
        if(edgeLowersDiff_.contains(i)){
            lower_ += edgeLowersDiff_[i];
        }
        else{
            lower_ += edgeLowersParent_->at(i);
        }

        if(edgeUppersDiff_.contains(i)){
            upper_ += edgeUppersDiff_[i];
        }
        else{
            upper_ += edgeUppersParent_->at(i);
        }
    }
}

void ArsGraphIntervalDiff::split(size_t idx, Base::Ptr intervLower, Base::Ptr intervUpper) {
    ARS_ASSERT_VAR1(idx != 0, idx);

    intervLower->setGraph(graph_);
    intervUpper->setGraph(graph_);

    double thetaLower = nodeLowersDiff_.contains(idx) ? 
        nodeLowersDiff_[idx] : nodeLowersParent_->at(idx);

    double thetaUpper = nodeUppersDiff_.contains(idx) ? 
        nodeUppersDiff_[idx] : nodeUppersParent_->at(idx);

    double thetaMid = 0.5 * (thetaLower + thetaUpper);

    intervLower->setNodeLowers(*nodeLowersParent_);
    intervLower->setNodeLowers(nodeLowersDiff_);
    intervLower->setNodeUppers(*nodeUppersParent_);
    intervLower->setNodeUppers(nodeUppersDiff_);
    intervLower->setNodeUpper(idx, thetaMid);

    intervUpper->setNodeLowers(*nodeLowersParent_);
    intervUpper->setNodeLowers(nodeLowersDiff_);
    intervUpper->setNodeUppers(*nodeUppersParent_);
    intervUpper->setNodeUppers(nodeUppersDiff_);
    intervUpper->setNodeLower(idx, thetaMid);

    intervLower->setEdgeLowers(*edgeLowersParent_);
    intervLower->setEdgeLowers(edgeLowersDiff_);
    intervLower->setEdgeUppers(*edgeUppersParent_);
    intervLower->setEdgeUppers(edgeUppersDiff_);

    intervUpper->setEdgeLowers(*edgeLowersParent_);
    intervUpper->setEdgeLowers(edgeLowersDiff_);
    intervUpper->setEdgeUppers(*edgeUppersParent_);
    intervUpper->setEdgeUppers(edgeUppersDiff_);

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

size_t ArsGraphIntervalDiff::size() const {
    size_t size = 0;
    size += sizeof(graph_);
    size += sizeof(nodeLowersParent_) + nodeLowersDiff_.size() *
        (sizeof(nodeLowersDiff_.begin()->first)+sizeof(nodeLowersDiff_.begin()->second));

    size += sizeof(nodeUppersParent_) + nodeUppersDiff_.size() *
        (sizeof(nodeUppersDiff_.begin()->first)+sizeof(nodeUppersDiff_.begin()->second));

    size += sizeof(edgeLowersParent_) + edgeLowersDiff_.size() *
        (sizeof(edgeLowersDiff_.begin()->first)+sizeof(edgeLowersDiff_.begin()->second));

    size += sizeof(edgeUppersParent_) + edgeUppersDiff_.size() *
        (sizeof(edgeUppersDiff_.begin()->first)+sizeof(edgeUppersDiff_.begin()->second));

    size += sizeof(lower_);
    size += sizeof(upper_);
    return size;
}

}; // namespace ars