#include <ars/BBTranslation.h>

namespace ars
{

    double distancePointBox(const Vector2 &p,
                            const Vector2 &boxMin,
                            const Vector2 &boxMax)
    {
        double dist = 0.0;
        double len;
        for (int d = 0; d < 2; ++d)
        {
            if (boxMin(d) <= p(d) && p(d) <= boxMax(d))
            {
                len = 0.0;
            }
            else if (p(d) < boxMin(d))
            {
                len = boxMin(d) - p(d);
            }
            else
            {
                len = p(d) - boxMax(d);
            }
            dist += len * len;
        }
        return dist;
    }

    Box::Box(const Vector2 &min, const Vector2 &max)
        : min_(min), max_(max), lower_(0.0), upper_(0.0) {}

    Box::Box(const Vector2 &min,
             const Vector2 &max,
             const VectorVector2 &ptsSrc,
             const VectorVector2 &ptsDst)
    {
        double dist, distMin, distUpper, distUpperMin;
        Vector2 boxMin, boxMax, boxMid;
        min_ = min;
        max_ = max;
        // computeBoundsNaive(ptsSrc, ptsDst);
        computeBoundsInlier(ptsSrc, ptsDst);
    }

    Box::~Box() {}

    void Box::computeBoundsNaive(const VectorVector2 &ptsSrc,
                                 const VectorVector2 &ptsDst)
    {
        double distLower, distLowerMin, distUpper, distUpperMin;
        Vector2 boxMin, boxMax, boxMid;
        lower_ = 0.0;
        upper_ = 0.0;
        for (int is = 0; is < ptsSrc.size(); ++is)
        {
            distLowerMin = 1e+6;
            distUpperMin = 1e+6;
            boxMin = ptsSrc[is] + min_;
            boxMax = ptsSrc[is] + max_;
            boxMid = ptsSrc[is] + 0.5 * (max_ + min_);
            for (int id = 0; id < ptsDst.size(); ++id)
            {
                distLower = distancePointBox(ptsDst[id], boxMin, boxMax);
                if (distLower < distLowerMin)
                {
                    distLowerMin = distLower;
                }
                distUpper = (ptsDst[id] - (ptsSrc[is] + boxMid)).squaredNorm();
                // ARS_VAR5(boxMid.transpose(), ptsSrc[is].transpose(),
                //          ptsDst[id].transpose(), distUpper, distUpperMin);
                if (distUpper < distUpperMin)
                {
                    distUpperMin = distUpper;
                }
            }
            lower_ += distLowerMin;
            upper_ += distUpperMin;
            // ARS_VAR4(distLowerMin, distUpperMin, lower_, upper_);
        }
    }

    void Box::computeBoundsInlier(const VectorVector2 &ptsSrc,
                                  const VectorVector2 &ptsDst)
    {
        const double eps = 0.1;
        Vector2 mid = 0.5 * (min_ + max_);
        Vector2 srcTransl;
        double dist, len;
        bool inlierFoundUpper, inlierFoundLower;

        len = 0.5 * (max_ - min_).maxCoeff(); // Half of Infinity norm
        lower_ = (double)ptsSrc.size();
        upper_ = (double)ptsSrc.size();
        ARS_VAR4(lower_, upper_, len, mid.transpose());
        for (int is = 0; is < ptsSrc.size(); ++is)
        {
            srcTransl = ptsSrc[is] + mid;
            inlierFoundLower = false;
            inlierFoundUpper = false;
            // ARS_VAR1(srcTransl.transpose());
            for (int id = 0; id < ptsDst.size() && !(inlierFoundLower && inlierFoundUpper); ++id)
            {
                // dist = (ptsDst[id] - srcTransl).norm();
                dist = (ptsDst[id] - srcTransl).cwiseAbs().maxCoeff(); // Infinity norm
                // ARS_VAR4(ptsDst[id].transpose(), dist, dist < eps, dist < eps + len);
                if (dist < eps)
                {
                    inlierFoundUpper = true;
                }
                if (dist < eps + len)
                {
                    inlierFoundLower = true;
                }
            }
            if (inlierFoundLower)
                lower_ -= 1.0;
            if (inlierFoundUpper)
                upper_ -= 1.0;
        }
    }

    std::ostream &operator<<(std::ostream &out, const Box &box)
    {
        out << "min [" << box.min_.transpose() << "] max [" << box.max_.transpose()
            << "] lower " << box.lower_ << " upper " << box.upper_;
        return out;
    }

    BBTranslation::BBTranslation() : res_(0.2) {}

    BBTranslation::~BBTranslation() {}

    void BBTranslation::compute()
    {
        Vector2 boxSplitMin, boxSplitMax;
        Vector2 translOpt;
        double scoreOpt, scoreTol;
        int iterNum;

        auto cmp = [](const Box &box1, const Box &box2)
        {
            return box1.lower_ > box2.lower_;
        };
        std::priority_queue<Box, std::vector<Box>, decltype(cmp)> prioqueue(cmp);

        scoreTol = 0.05; // TODO: allow setting the value of scoreTol
        Box boxCur(translMin_, translMax_, ptsSrc_, ptsDst_);
        prioqueue.push(boxCur);
        scoreOpt = prioqueue.top().upper_;
        translOpt = 0.5 * (boxCur.min_ + boxCur.max_);
        ARS_VAR2(boxCur, scoreOpt);
        iterNum = 0;
        while (!prioqueue.empty() && iterNum < 100000)
        {
            boxCur = prioqueue.top();
            prioqueue.pop();

            std::cout << "\n---\niteration " << iterNum << " queue size "
                      << prioqueue.size() << std::endl;
            ARS_PRINT("boxCur " << boxCur << ", score optimum " << scoreOpt);
            // ARS_VAR4(boxCur.upper_, boxCur.lower_, scoreTol * scoreOpt,
            //          boxCur.upper_ - boxCur.lower_ <= scoreTol * scoreOpt);
            if (scoreOpt - boxCur.lower_ <= scoreTol * scoreOpt)
            {
                ARS_PRINT("STOP");
                break;
            }

            // Splits the current box into 2^DIM parts
            if ((boxCur.max_ - boxCur.min_).maxCoeff() > res_)
            {
                for (int j = 0; j < SPLIT_NUM; ++j)
                {
                    for (int d = 0; d < DIM; ++d)
                    {
                        // ARS_VAR4(j, d, (1 << d), j & (1 << d));
                        if (j & (1 << d))
                        {
                            boxSplitMin(d) = 0.5 * (boxCur.min_(d) + boxCur.max_(d));
                            boxSplitMax(d) = boxCur.max_(d);
                        }
                        else
                        {
                            boxSplitMin(d) = boxCur.min_(d);
                            boxSplitMax(d) = 0.5 * (boxCur.min_(d) + boxCur.max_(d));
                        }
                        // ARS_VAR2(boxSplitMin(d), boxSplitMax(d));
                    }
                    Box boxNew(boxSplitMin, boxSplitMax, ptsSrc_, ptsDst_);
                    ARS_VAR1(boxNew);

                    if (boxNew.upper_ < scoreOpt)
                    {
                        scoreOpt = boxNew.upper_;
                        translOpt = 0.5 * (boxNew.min_ + boxNew.max_);
                        ARS_PRINT("UPDATE optimum " << scoreOpt << " in "
                                                    << translOpt.transpose());
                    }

                    if (boxNew.lower_ < scoreOpt)
                    {
                        prioqueue.push(boxNew);
                    }
                }
            }
            iterNum++;
        }
        ARS_PRINT("OPTIMUM " << scoreOpt << " in " << translOpt.transpose());
    }

    void BBTranslation::setTranslMinMax(const ars::Vector2 &translMin,
                                        const ars::Vector2 &translMax)
    {
        translMin_ = translMin;
        translMax_ = translMax;
        ARS_VAR2(translMin_.transpose(), translMax_.transpose())
    }

    void BBTranslation::setResolution(const double r)
    {
        res_ = r;
    }

    void BBTranslation::setPtsSrc(const ars::VectorVector2 &pts)
    {
        ptsSrc_ = pts;
    }

    void BBTranslation::setPtsDst(const ars::VectorVector2 &pts)
    {
        ptsDst_ = pts;
    }

    void BBTranslation::setPts(const ars::VectorVector2 &ptsS,
                               const ars::VectorVector2 &ptsD)
    {
        ptsSrc_ = ptsS;
        ptsDst_ = ptsD;
    }

} // namespace ars