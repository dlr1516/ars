#include <ars/ArsGraphSolver.h>

namespace ars {

ArsGraphSolver::ArsGraphSolver() : graph_(nullptr), lower_(0.0), upper_(0.0) {}

ArsGraphSolver::ArsGraphSolver(ArsGraphPtr& graph)
    : graph_(graph), lower_(0.0), upper_(0.0) {}

ArsGraphSolver::~ArsGraphSolver() {}

void ArsGraphSolver::setGraph(ArsGraphPtr& graph) {
    graph_ = graph;
}

bool ArsGraphSolver::solve(std::vector<double>& solution, double& cost) {
    LeastUpperBoundFirstQueue queue;
    ArsGraphIntervalFull::Ptr initial =
        std::make_shared<ArsGraphIntervalFull>(graph_);

    initial->init(graph_);

    queue.push(initial);
    lower_ = initial->getLower();
    upper_ = initial->getUpper();

    while (!queue.empty()) {
        ArsGraphIntervalPtr curr = queue.top();
        queue.pop();
        if (curr->getUpper() >= lower_) {
            if (lower_ < curr->getLower()) {
                lower_ = curr->getLower();
                upper_ = curr->getUpper();
                solution_.angleLower = curr->getNodesLower();
                solution_.angleUpper = curr->getNodesUpper();
            }
            // if (!curr->isAtomic()) {
            //     ArsGraphIntervalPtr left = curr->splitLeft();
            //     ArsGraphIntervalPtr right = curr->splitRight();
            //     if (left->isValid()) {
            //         queue.push(left);
            //     }
            //     if (right->isValid()) {
            //         queue.push(right);
            //     }
            // }
        }
    }
    cost = lower_;
    return true;
}

}  // namespace ars