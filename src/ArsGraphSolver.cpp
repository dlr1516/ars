#include <ars/ArsGraphSolver.h>

namespace ars {

ArsGraphSolver::ArsGraphSolver() : graph_(nullptr), lower_(0.0), upper_(0.0), xtol_(M_PI / 180.0) {}

ArsGraphSolver::ArsGraphSolver(ArsGraphPtr& graph)
    : graph_(graph), lower_(0.0), upper_(0.0), xtol_(M_PI / 180.0) {}

ArsGraphSolver::ArsGraphSolver(ArsGraphPtr& graph, double xtol)
    : graph_(graph), lower_(0.0), upper_(0.0), xtol_(xtol) {}

ArsGraphSolver::~ArsGraphSolver() {}

void ArsGraphSolver::setGraph(ArsGraphPtr& graph) {
    graph_ = graph;
}

void ArsGraphSolver::setXTol(double xtol) {
    xtol_ = xtol;
}

bool ArsGraphSolver::solve(std::vector<double>& solution, double& cost, bool useDiff) {
    LeastUpperBoundFirstQueuePtr queue(new LeastUpperBoundFirstQueue);
    ArsGraphIntervalFull::Ptr initial(new ArsGraphIntervalFull(graph_));
    initial->init(graph_);
    queue->push(initial);

    solution.clear();

    lower_ = initial->getLowerBound();
    upper_ = initial->getUpperBound();
    solution_.anglesLower = initial->getNodeLowers();
    solution_.anglesUpper = initial->getNodeUppers();

    while (!queue->empty()) {
        ArsGraphIntervalPtr curr = queue->top();
        queue->pop();
        if (curr->getUpperBound() >= lower_) {
            if (lower_ < curr->getLowerBound()) {
                lower_ = curr->getLowerBound();
                upper_ = curr->getUpperBound();
                solution_.anglesLower = curr->getNodeLowers();
                solution_.anglesUpper = curr->getNodeUppers();
            }
            std::vector<NodeInterval> validNodes;
            if(!checkInterval(curr, validNodes)){
                ars::ArsGraphInterval::Ptr intervLower;
                ars::ArsGraphInterval::Ptr intervUpper;
                if (useDiff && curr->getDiffSize() <= curr->getEdgeNum()/3) {
                    intervLower.reset(new ars::ArsGraphIntervalDiff);
                    intervUpper.reset(new ars::ArsGraphIntervalDiff);
                }
                else {
                    intervLower.reset(new ars::ArsGraphIntervalFull);
                    intervUpper.reset(new ars::ArsGraphIntervalFull);
                }

                int node = validNodes[0].idx;
                for(int i = 1; i < validNodes.size(); i++){
                    if (validNodes[node].xWidth < validNodes[i].xWidth){
                        node = validNodes[i].idx;
                        break;
                    }
                }

                curr->split(node, intervLower, intervUpper);
                queue->push(intervLower);
                queue->push(intervUpper);
            }
        }
    }
    solution.resize(solution_.anglesLower.size());
    for(int i = 0; i < solution_.anglesLower.size(); i++){
        double angleLower = solution_.anglesLower[i];
        double angleUpper = solution_.anglesUpper[i];
        solution[i] = (angleLower + angleUpper)/2;
    }
    cost = lower_;
    return true;
}

bool ArsGraphSolver::solve(std::vector<double>& solution, double& cost, 
    Statistics& stats, bool useDiff) {

    LeastUpperBoundFirstQueuePtr queue(new LeastUpperBoundFirstQueue);
    ArsGraphIntervalFull::Ptr initial(new ArsGraphIntervalFull(graph_));
    initial->init(graph_);
    queue->push(initial);

    solution.clear();
    stats.createdNodes = 1;
    stats.minIntervalSize = initial->size();
    stats.maxIntervalSize = initial->size();
    stats.avgIntervalSize = initial->size();

    lower_ = initial->getLowerBound();
    upper_ = initial->getUpperBound();
    solution_.anglesLower = initial->getNodeLowers();
    solution_.anglesUpper = initial->getNodeUppers();

    while (!queue->empty()) {
        ArsGraphIntervalPtr curr = queue->top();
        queue->pop();
        if (curr->getUpperBound() >= lower_) {
            if (lower_ < curr->getLowerBound()) {
                lower_ = curr->getLowerBound();
                upper_ = curr->getUpperBound();
                solution_.anglesLower = curr->getNodeLowers();
                solution_.anglesUpper = curr->getNodeUppers();
            }
            std::vector<NodeInterval> validNodes;
            if(!checkInterval(curr, validNodes)){
                ars::ArsGraphInterval::Ptr intervLower;
                ars::ArsGraphInterval::Ptr intervUpper;
                if (useDiff && curr->getDiffSize() <= curr->getEdgeNum()/3) {
                    intervLower.reset(new ars::ArsGraphIntervalDiff);
                    intervUpper.reset(new ars::ArsGraphIntervalDiff);
                }
                else {
                    intervLower.reset(new ars::ArsGraphIntervalFull);
                    intervUpper.reset(new ars::ArsGraphIntervalFull);
                }

                int idx = 0;
                for(int i = 1; i < validNodes.size(); i++){
                    if (validNodes[idx].xWidth < validNodes[i].xWidth){
                        idx = i;
                        break;
                    }
                }
                int node = validNodes[idx].idx;

                curr->split(node, intervLower, intervUpper);
                queue->push(intervLower);
                queue->push(intervUpper);
                stats.createdNodes += 2;
                size_t lowerSize = intervLower->size();
                size_t upperSize = intervUpper->size();

                if(lowerSize < stats.minIntervalSize){
                    stats.minIntervalSize = lowerSize;
                }
                if(lowerSize > stats.maxIntervalSize){
                    stats.maxIntervalSize = lowerSize;
                }

                if(upperSize < stats.minIntervalSize){
                    stats.minIntervalSize = upperSize;
                }
                if(upperSize > stats.maxIntervalSize){
                    stats.maxIntervalSize = upperSize;
                }

                stats.avgIntervalSize = 
                    (stats.avgIntervalSize * (stats.createdNodes-2) 
                        + lowerSize + upperSize) / stats.createdNodes;
            }
        }
    }
    solution.resize(solution_.anglesLower.size());
    for(int i = 0; i < solution_.anglesLower.size(); i++){
        double angleLower = solution_.anglesLower[i];
        double angleUpper = solution_.anglesUpper[i];
        solution[i] = (angleLower + angleUpper)/2;
    }
    cost = lower_;
    return true;
}

bool ArsGraphSolver::checkInterval(ArsGraphIntervalPtr interval, std::vector<NodeInterval>& validNodes){
    validNodes.clear();

    for(size_t i = 1; i < interval->getNodeNum(); i++){
        double xLower = interval->nodeLower(i);
        double xUpper = interval->nodeUpper(i);
        double xWidth = xUpper - xLower;
        if(xWidth > xtol_){
            validNodes.push_back(NodeInterval(i, xLower, xUpper, xWidth));
        }
    }
    return validNodes.empty();
}

}  // namespace ars