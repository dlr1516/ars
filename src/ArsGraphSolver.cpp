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

void ArsGraphSolver::setXTol(double xtol){
    xtol_ = xtol;
}

bool ArsGraphSolver::solve(std::vector<double>& solution, double& cost) {
    LeastUpperBoundFirstQueue queue;
    ArsGraphIntervalFull::Ptr initial =
        std::make_shared<ArsGraphIntervalFull>(graph_);

    initial->init(graph_);

    queue.push(initial);
    lower_ = initial->getLowerBound();
    upper_ = initial->getUpperBound();

    while (!queue.empty()) {
        ArsGraphIntervalPtr curr = queue.top();
        queue.pop();
        if (curr->getUpperBound() >= lower_) {
            if (lower_ < curr->getLowerBound()) {
                lower_ = curr->getLowerBound();
                upper_ = curr->getUpperBound();
                solution_.anglesLower = curr->getNodesLower();
                solution_.anglesUpper = curr->getNodesUpper();
            }
            std::vector<NodeInterval> validNodes;
            if(!checkInterval(curr, validNodes)){
                ars::ArsGraphIntervalFull::Ptr intervLower;
                ars::ArsGraphIntervalFull::Ptr intervUpper;

                int node = validNodes[0].idx;
                for(int i = 1; i < validNodes.size(); i++){
                    if (validNodes[node].xWidth < validNodes[i].xWidth){
                        node = validNodes[i].idx;
                        break;
                    }
                }

                curr->split(node, intervLower, intervUpper);
                queue.push(intervLower);
                queue.push(intervUpper);
            }
        }
    }
    cost = lower_;
    return true;
}

bool ArsGraphSolver::checkInterval(ArsGraphIntervalPtr interval, std::vector<NodeInterval>& validNodes){
    validNodes.clear();

    for(int i = 1; i < interval->getNodeNum(); i++){
        double xLower = interval->nodeLower(i);
        double xUpper = interval->nodeUpper(i);
        double xWidth = xUpper - xLower;
        if(xWidth > xtol_){
            validNodes.push_back(NodeInterval(xLower, xUpper, xWidth));
        }
    }
    return validNodes.empty();
}

}  // namespace ars