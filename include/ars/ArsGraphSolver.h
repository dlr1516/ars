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
            return (ib0->getUpper() < ib1->getUpper());
        }
    };

    using LeastUpperBoundFirstQueue =
        std::priority_queue<ArsGraphIntervalPtr,
                            std::vector<ArsGraphIntervalPtr>,
                            UpperBoundLess>;

    struct Solution {
        std::vector<double> angleLower;
        std::vector<double> angleUpper;
    };

    ArsGraphSolver();

    ArsGraphSolver(ArsGraphPtr& graph);

    virtual ~ArsGraphSolver();

    void setGraph(ArsGraphPtr& graph);

    bool solve(std::vector<double>& solution, double& cost);

   protected:
    ArsGraphPtr graph_;
    double lower_;
    double upper_;
    Solution solution_;
};

}  // namespace ars

#endif  // ARS_GRAPH_SOLVER_H