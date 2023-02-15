#include <ars/BBTranslation.h>

namespace ars {

double distancePointBox(const Vector2& p,
                        const Vector2& boxMin,
                        const Vector2& boxMax) {
    double dist = 0.0;
    double len;
    for (int d = 0; d < 2; ++d) {
        if (boxMin(d) <= p(d) && p(d) <= boxMax(d)) {
            len = 0.0;
        } else if (p(d) < boxMin(d)) {
            len = boxMin(d) - p(d);
        } else {
            len = p(d) - boxMax(d);
        }
        dist += len * len;
    }
    return dist;
}

Box::Box(const Vector2& min, const Vector2& max)
    : min_(min), max_(max), lower_(0.0), upper_(0.0) {}

Box::Box(const Vector2& min,
         const Vector2& max,
         const VectorVector2& ptsSrc,
         const VectorVector2& ptsDst) {
    double dist, distMin, distUpper, distUpperMin;
    Vector2 boxMin, boxMax, boxMid;
    min_ = min;
    max_ = max;
    computeBounds(ptsSrc, ptsDst);
}

Box::~Box() {}

void Box::computeBounds(const VectorVector2& ptsSrc,
                        const VectorVector2& ptsDst) {
    double distLower, distLowerMin, distUpper, distUpperMin;
    Vector2 boxMin, boxMax, boxMid;
    lower_ = 0.0;
    upper_ = 0.0;
    for (int is = 0; is < ptsSrc.size(); ++is) {
        distLowerMin = 1e+6;
        distUpperMin = 1e+6;
        boxMin = ptsSrc[is] + min_;
        boxMax = ptsSrc[is] + max_;
        boxMid = ptsSrc[is] + 0.5 * (max_ + min_);
        for (int id = 0; id < ptsDst.size(); ++id) {
            distLower = distancePointBox(ptsDst[id], boxMin, boxMax);
            if (distLower < distLowerMin) {
                distLowerMin = distLower;
            }
            distUpper = (ptsDst[id] - boxMid).squaredNorm();
            if (distUpper < distUpperMin) {
                distUpperMin = distUpper;
            }
        }
        lower_ += distLowerMin;
        upper_ += distUpperMin;
    }
}

BBTranslation::BBTranslation() {}

BBTranslation::~BBTranslation() {}

void BBTranslation::compute() {
    Vector2 boxSplitMin, boxSplitMax;
    Vector2 translOpt;
    double scoreOpt, scoreTol;

    auto cmp = [](const Box& box1, const Box& box2) {
        return box1.upper_ < box2.upper_;
    };
    std::priority_queue<Box, std::vector<Box>, decltype(cmp)> prioqueue(cmp);

    scoreTol = 0.05;  // TODO: allow setting the value of scoreTol
    prioqueue.push(Box(translMin_, translMax_, ptsSrc_, ptsDst_));
    scoreOpt = prioqueue.top().upper_;
    while (!prioqueue.empty()) {
        Box boxCur = prioqueue.top();
        prioqueue.pop();

        if (boxCur.upper_ - scoreOpt <= scoreTol * scoreOpt) {
            break;
        }

        if (boxCur.upper_ < scoreOpt) {
            scoreOpt = boxCur.upper_;
            translOpt = 0.5 * (boxCur.min_ + boxCur.max_);
        }

        // Splits the current box into 2^DIM parts
        for (int j = 0; j < (1 >> 2); ++j) {
            for (int d = 0; d < 2; ++d) {
                if (j & (1 << d)) {
                    boxSplitMin(d) = 0.5 * (boxCur.min_(d) + boxCur.max_(d));
                    boxSplitMax(d) = boxCur.max_(d);
                } else {
                    boxSplitMin(d) = boxCur.min_(d);
                    boxSplitMax(d) = 0.5 * (boxCur.min_(d) + boxCur.max_(d));
                }
            }
            Box boxNew(boxSplitMin, boxSplitMax, ptsSrc_, ptsDst_);
            if (boxNew.lower_ < scoreOpt) {
                prioqueue.push(boxNew);
            }
        }
    }
}

void BBTranslation::setTranslMinMax(const ars::Vector2& translMin,
                                    const ars::Vector2& translMax) {
    translMin_ = translMin;
    translMax_ = translMax;
}

void BBTranslation::setPtsSrc(const ars::VectorVector2& pts) {
    ptsSrc_ = pts;
}

void BBTranslation::setPtsDst(const ars::VectorVector2& pts) {
    ptsDst_ = pts;
}

void BBTranslation::setPts(const ars::VectorVector2& ptsS,
                           const ars::VectorVector2& ptsD) {
    ptsSrc_ = ptsS;
    ptsDst_ = ptsD;
}

}  // namespace ars