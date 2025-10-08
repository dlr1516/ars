#ifndef ARS_GRAPH_SOLVER_H
#define ARS_GRAPH_SOLVER_H

#include <queue>

#include <ars/ArsGraph.h>
#include <ars/ArsGraphInterval.h>

namespace ars {

class ArsGraphSolver {
   public:
    using Self = ArsGraphSolver;
    using Ptr = std::shared_ptr<Self>;
    using ArsGraphPtr = ArsGraph::Ptr;
    using ArsGraphIntervalPtr = ArsGraphInterval::Ptr;

    /**
     * Comparator of intervals to order interval with smaller upper bound
     * before.
     */
    struct UpperBoundLess {
        bool operator()(ArsGraphIntervalPtr& ib0,
                        ArsGraphIntervalPtr& ib1) const {
            return (ib0->getUpperBound() < ib1->getUpperBound());
        }
    };

    struct NodeInterval {
        size_t idx;
        double xLower;
        double xUpper;
        double xWidth;

        NodeInterval(size_t idx, double xL, double xU, double xW) 
            : idx(idx), xLower(xL), xUpper(xU), xWidth(xW) {}
    };

    using LeastUpperBoundFirstQueue =
        std::priority_queue<ArsGraphIntervalPtr,
                            std::vector<ArsGraphIntervalPtr>,
                            UpperBoundLess>;

    struct Solution {
        std::vector<double> anglesLower;
        std::vector<double> anglesUpper;
    };

    ArsGraphSolver();

    ArsGraphSolver(ArsGraphPtr& graph);
    ArsGraphSolver(ArsGraphPtr& graph, double xtol);

    virtual ~ArsGraphSolver();

    void setGraph(ArsGraphPtr& graph);
    void setXTol(double xtol);

    bool solve(std::vector<double>& solution, double& cost);

   protected:
    ArsGraphPtr graph_;
    double lower_;
    double upper_;
    Solution solution_;
    double xtol_;

    bool checkInterval(ArsGraphIntervalPtr interval, std::vector<NodeInterval>& validIndices);
};

}  // namespace ars

#endif  // ARS_GRAPH_SOLVER_H