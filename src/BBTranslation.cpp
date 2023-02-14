#include <ars/BBTranslation.h>

namespace ars
{

    double distancePointBox(const Vector2& p, const Vector2 &boxMin, const Vector2 &boxMax) {
        double dist = 0.0;
        double len;
        for (int d = 0; d < 2; ++d) {
            if (boxMin(d) <= p(d) && p(d) <= boxMax(d)) {
                len = 0.0;
            }
            else if (p(d) < boxMin(d)) {
                len = boxMin(d) - p(d); 
            }
            else {
                len = p(d) - boxMax(d);
            }
            dist += len * len;
        }
        return dist;
    }

    BBTranslation::BBTranslation();

    BBTranslation::~BBTranslation();

    void BBTranslation::compute() {
        auto cmp = [](const Box& box1,const Box& box2) {
            return box1.lower < box2.lower;
        };
        std::priority_queue<Box, std::vector<Box>, decltype(cmp)> prioqueue; 


    }

    void BBTranslation::setTranslMinMax(const ars::Vector2 &translMin, const ars::Vector2 &translMax) {
        translMin_ = translMin;
        translMax_ = translMax;
    }

    void BBTranslation::setPtsSrc(const ars::VectorVector2 &pts)
    {
        ptsSrc_ = pts;
    }

    void BBTranslation::setPtsDst(const ars::VectorVector2 &pts)
    {
        ptsDst_ = pts;
    }

    void BBTranslation::setPts(const ars::VectorVector2 &ptsS, const ars::VectorVector2 &ptsD)
    {
        ptsSrc_ = ptsS;
        ptsDst_ = ptsD;
    }

}