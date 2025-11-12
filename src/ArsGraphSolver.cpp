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

bool ArsGraphSolver::initialSolutionFromInitialAndTree(const ArsGraphInterval::Ptr& initial) {
    lower_ = initial->getLowerBound();
    upper_ = initial->getUpperBound();
    solution_.anglesLower = initial->getNodeLowers();
    solution_.anglesUpper = initial->getNodeUppers();

    std::vector<int> added;
    added.push_back(0);

    std::queue<int> explorable;
    explorable.push(0);

    while(!explorable.empty()){
        int isrc = explorable.front();
        auto node = graph_->nodes()[isrc];
        explorable.pop();
        for(int eidx : node.incidents) {
            auto edge = graph_->edges()[eidx];
            int mult = 1;
            int idst;

            if (isrc == edge.isrc){      
                idst = edge.idst;
            }
            else{
                idst = edge.isrc;
                mult = -1;
            }

            if(std::find(added.begin(), added.end(), idst) != added.end())
                continue;

            double tMax, fLow, fUp;
            FourierOptimizerBB1D fopt(edge.coeffs);

            fopt.setXTolerance(xtol_);
            fopt.enableXTolerance(true);
            fopt.enableYTolerance(false);

            fopt.findGlobalMax(.0, M_PI, tMax, fLow, fUp);

            lower_ += (fLow - initial->edgeLower(eidx));
            upper_ += (fUp - initial->edgeUpper(eidx));
            double sol = solution_.anglesLower[isrc] + (mult*tMax);
            solution_.anglesLower[idst] = sol;
            solution_.anglesUpper[idst] = sol;

            explorable.push(idst);
            added.push_back(idst);
        }   
    }
    return added.size() == graph_->nodes().size();
}

bool ArsGraphSolver::initialSolutionFromTree() {
    lower_ = .0;
    upper_ = .0;

    solution_.anglesLower.resize(graph_->getNodeNum());
    solution_.anglesUpper.resize(graph_->getNodeNum());
    solution_.anglesLower[0] = .0;
    solution_.anglesUpper[0] = .0;

    std::vector<int> added;
    added.push_back(0);

    std::queue<int> explorable;
    explorable.push(0);

    while(!explorable.empty()){
        int isrc = explorable.front();
        auto node = graph_->nodes()[isrc];
        explorable.pop();
        for(int eidx : node.incidents) {
            auto edge = graph_->edges()[eidx];
            int mult = 1;
            int idst;

            if (isrc == edge.isrc){      
                idst = edge.idst;
            }
            else{
                idst = edge.isrc;
                mult = -1;
            }

            if(std::find(added.begin(), added.end(), idst) != added.end()){
                if (isrc == edge.isrc){
                    double fLow, fUp;

                    double t = solution_.anglesLower[idst] - solution_.anglesLower[isrc];

                    findLUFourier(edge.coeffs, t - xtol_*.5, t + xtol_*.5, fLow, fUp);

                    lower_ += fLow;
                    upper_ += fUp;
                }
                continue;
            }

            double tMax, fLow, fUp;
            FourierOptimizerBB1D fopt(edge.coeffs);

            fopt.setXTolerance(xtol_);
            fopt.enableXTolerance(true);
            fopt.enableYTolerance(false);

            fopt.findGlobalMax(.0, M_PI, tMax, fLow, fUp);

            lower_ += fLow;
            upper_ += fUp;
            double sol = solution_.anglesLower[isrc] + (mult*tMax);
            solution_.anglesLower[idst] = sol;
            solution_.anglesUpper[idst] = sol;

            explorable.push(idst);
            added.push_back(idst);
        }   
    }
    return added.size() == graph_->nodes().size();
}

bool ArsGraphSolver::initialSolutionFromGraph() {
    lower_ = 0;
    upper_ = 0;

    solution_.anglesLower.resize(graph_->getNodeNum());
    solution_.anglesUpper.resize(graph_->getNodeNum());
    solution_.anglesLower[0] = .0;
    solution_.anglesUpper[0] = .0;

    for(int i = 0; i < graph_->getNodeNum(); i++){
        auto node = graph_->nodes()[i];
        for(int eidx : node.incidents) {
            auto edge = graph_->edges()[eidx];
            int idst = edge.idst;

            if (i != edge.isrc){      
                continue;
            }

            double tMax, fLow, fUp;
            FourierOptimizerBB1D fopt(edge.coeffs);

            fopt.setXTolerance(xtol_);
            fopt.enableXTolerance(true);
            fopt.enableYTolerance(false);

            fopt.findGlobalMax(.0, M_PI, tMax, fLow, fUp);

            lower_ += fLow;
            upper_ += fUp;
            double sol = solution_.anglesLower[i] + (tMax);
            solution_.anglesLower[idst] = sol;
            solution_.anglesUpper[idst] = sol;
        }
        
    }
    return true;
}

bool ArsGraphSolver::solve(std::vector<double>& solution, double& cost, bool useDiff) {    
    LeastUpperBoundFirstQueuePtr queue(new LeastUpperBoundFirstQueue);
    ArsGraphIntervalFull::Ptr initial(new ArsGraphIntervalFull(graph_));
    initial->init(graph_);
    queue->push(initial);

    if(solution_.anglesLower.size() != graph_->getNodeNum()){
        lower_ = initial->getLowerBound();
        upper_ = initial->getUpperBound();
        solution_.anglesLower = initial->getNodeLowers();
        solution_.anglesUpper = initial->getNodeUppers();
    }
    std::cout << "Initial solution: " << std::endl;
    for(int i = 0; i < solution_.anglesLower.size(); i++){
        std::cout << "Node " << i << ": " 
            << solution_.anglesLower[i]*(180.0/M_PI) <<
            " - " << solution_.anglesUpper[i]*(180.0/M_PI) << std::endl; 
    }

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

    stats.createdNodes = 1;
    stats.minIntervalSize = initial->size();
    stats.maxIntervalSize = initial->size();
    stats.avgIntervalSize = initial->size();

    bool updated = false;

    if(solution_.anglesLower.size() != graph_->getNodeNum()){
        lower_ = initial->getLowerBound();
        upper_ = initial->getUpperBound();
        solution_.anglesLower = initial->getNodeLowers();
        solution_.anglesUpper = initial->getNodeUppers();
    }
    std::cout << "Initial solution: " << std::endl;
    for(int i = 0; i < solution_.anglesLower.size(); i++){
        std::cout << "Node " << i << ": " 
            << solution_.anglesLower[i]*(180.0/M_PI) << 
            " - " << solution_.anglesUpper[i]*(180.0/M_PI) << std::endl; 
    }

    while (!queue->empty()) {
        ArsGraphIntervalPtr curr = queue->top();
        queue->pop();
        if (curr->getUpperBound() >= lower_) {
            if (lower_ < curr->getLowerBound()) {
                lower_ = curr->getLowerBound();
                upper_ = curr->getUpperBound();
                solution_.anglesLower = curr->getNodeLowers();
                solution_.anglesUpper = curr->getNodeUppers();
                if(!updated){
                    std::cout << "Updated sol after: " << stats.createdNodes <<std::endl;
                    updated = true;
                }
            }
            std::vector<NodeInterval> validNodes;
            if(!checkInterval(curr, validNodes)){
                ars::ArsGraphInterval::Ptr intervLower;
                ars::ArsGraphInterval::Ptr intervUpper;
                if (useDiff && curr->getDiffSize() <= curr->getEdgeNum()*0.5) {
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

bool ArsGraphSolver::solveWithStationary(std::vector<double>& solution, double& cost, 
    Statistics& stats, bool useDiff) {
    
    LeastUpperBoundFirstQueuePtr queue(new LeastUpperBoundFirstQueue);
    ArsGraphIntervalFull::Ptr initial(new ArsGraphIntervalFull(graph_));
    initial->initWithStationary(graph_);
    queue->push(initial);

    stats.createdNodes = 1;
    stats.minIntervalSize = initial->size();
    stats.maxIntervalSize = initial->size();
    stats.avgIntervalSize = initial->size();

    bool updated = false;

    if(solution_.anglesLower.size() != graph_->getNodeNum()){
        lower_ = initial->getLowerBound();
        upper_ = initial->getUpperBound();
        solution_.anglesLower = initial->getNodeLowers();
        solution_.anglesUpper = initial->getNodeUppers();
    }
    std::cout << "Initial solution: " << std::endl;
    for(int i = 0; i < solution_.anglesLower.size(); i++){
        std::cout << "Node " << i << ": " 
            << solution_.anglesLower[i]*(180.0/M_PI) << 
            " - " << solution_.anglesUpper[i]*(180.0/M_PI) << std::endl; 
    }

    while (!queue->empty()) {
        ArsGraphIntervalPtr curr = queue->top();
        queue->pop();
        if (curr->getUpperBound() >= lower_) {
            if (lower_ < curr->getLowerBound()) {
                lower_ = curr->getLowerBound();
                upper_ = curr->getUpperBound();
                solution_.anglesLower = curr->getNodeLowers();
                solution_.anglesUpper = curr->getNodeUppers();
                if(!updated){
                    std::cout << "Updated sol after: " << stats.createdNodes <<std::endl;
                    updated = true;
                }
            }
            std::vector<NodeInterval> validNodes;
            if(!checkInterval(curr, validNodes)){
                ars::ArsGraphInterval::Ptr intervLower;
                ars::ArsGraphInterval::Ptr intervUpper;
                if (useDiff && curr->getDiffSize() <= curr->getEdgeNum()*0.5) {
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

                curr->splitWithStationary(node, intervLower, intervUpper);
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
} // namespace ars