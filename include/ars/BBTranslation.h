#ifndef ARS_BBTRANSLATION_H_
#define ARS_BBTRANSLATION_H_

#include <ars/definitions.h>

#include <queue>

#include <Eigen/Dense>
#include <Eigen/Core>

namespace ars
{

    inline double distancePointBox(const Vector2 &p,
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

    struct Box
    {
        Vector2 min_;
        Vector2 max_;
        double lower_;
        double upper_;
        double eps_;

        Box(const Vector2 &min, const Vector2 &max, const double eps)
            : min_(min), max_(max), lower_(0.0), upper_(0.0), eps_(eps) {}

        Box(const Vector2 &min,
                 const Vector2 &max,
                 const VectorVector2 &ptsSrc,
                 const VectorVector2 &ptsDst,
                 const double eps)
        {
            double dist, distMin, distUpper, distUpperMin;
            Vector2 boxMin, boxMax, boxMid;
            min_ = min;
            max_ = max;
            eps_ = eps;
            // computeBoundsNaive(ptsSrc, ptsDst);
            computeBoundsInlier(ptsSrc, ptsDst);
        }

        ~Box() {}

        void computeBoundsNaive(const VectorVector2 &ptsSrc,
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

        void computeBoundsInlier(const VectorVector2 &ptsSrc,
                                      const VectorVector2 &ptsDst)
        {
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
                    // ARS_VAR4(ptsDst[id].transpose(), dist, dist < eps_, dist < eps_ + len);
                    if (dist < eps_)
                    {
                        inlierFoundUpper = true;
                    }
                    if (dist < eps_ + len)
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
    };

    inline std::ostream &operator<<(std::ostream &out, const Box &box)
    {
        out << "min [" << box.min_.transpose() << "] max [" << box.max_.transpose()
            << "] lower " << box.lower_ << " upper " << box.upper_;
        return out;
    }

    class BBTranslation
    {
    public:
        static constexpr int DIM = 2;
        static constexpr int SPLIT_NUM = (1 << DIM);

        /**
         * @brief Default constructor for a new BBTranslation object
         */
        BBTranslation();

        /**
         * @brief Default destructor for BBTranslation objects
         */
        virtual ~BBTranslation();

        /**
         * @brief Main method
         */
        void compute(Vector2 &translOpt);

        /**
         * @brief Set the interval search of translation
         */
        void setTranslMinMax(const ars::Vector2 &translMin,
                             const ars::Vector2 &translMax);

        /**
         * @brief Set minimum box size
         */
        void setResolution(const double r);

        /**
         * @brief Set points src
         */
        void setPtsSrc(const ars::VectorVector2 &pts);

        /**
         * @brief Set points dst
         */
        void setPtsDst(const ars::VectorVector2 &pts);

        /**
         * @brief Set pts src and dst
         */
        void setPts(const ars::VectorVector2 &ptsS, const ars::VectorVector2 &ptsD);

        /**
         * @brief Set epsilon param for lower bound computation
         */
        void setEps(const double eps);

        /**
         * @brief Set max number of iteration before B&B alg stops
         */
        void setNumMaxIterations(const int nmi);

    private:
        ars::Vector2 translMin_;
        ars::Vector2 translMax_;

        ars::VectorVector2 ptsSrc_;
        ars::VectorVector2 ptsDst_;

        double res_;
        double eps_;

        int numMaxIter_;
    };
} // namespace ars

#endif /*ARS_BBTRANSLATION_H_*/
