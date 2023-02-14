#ifndef ARS_BBTRANSLATION_H_
#define ARS_BBTRANSLATION_H_

#include <ars/definitions.h>

#include <queue>

#include <Eigen/Dense>

namespace ars
{

    double distancePointBox(const Vector2& p, const Vector2 &boxMin, const Vector2 &boxMax);

    struct Box
    {
        Vector2 min_;
        Vector2 max_;
        double lower_;
        double upper_;

        Box(const Vector2 &min, const Vector2 &max, const VectorVector2 &ptsSrc, const VectorVector2 &ptsDst)
        {
            double dist, distMin, distUpper, distUpperMin;
            Vector2 boxMin, boxMax, boxMid;
            min_ = min;
            max_ = max;
            lower_ = 0.0;
            upper_ = 0.0;
            for (int is = 0; is < ptsSrc.size(); ++is)
            {
                distMin = 1e+6;
                distUpperMin = 1e+6;
                boxMin = ptsSrc[is] + min_; 
                boxMax = ptsSrc[is] + max_;
                boxMid = ptsSrc[is] + 0.5 * (max_ + min_); 
                for (int id = 0; id < ptsDst.size(); ++id)
                {
                    dist = distancePointBox(ptsDst[id],boxMin,boxMax);
                    if (dist < distMin) {
                        distMin = dist;
                    }
                    distUpper = (ptsDst[id] - boxMid).squaredNorm();
                    if (distUpper < distUpperMin) {
                        distUpperMin = distUpper;
                    }
                }
                lower_ += distMin;
                upper_ += distUpperMin; 
            }
        }
    };

    class BBTranslation
    {
    public:
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
        void compute();

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

    private:
        ars::Vector2 translMin_;
        ars::Vector2 translMax_;

        ars::VectorVector2 ptsSrc_;
        ars::VectorVector2 ptsDst_;
    };
}

#endif /*ARS_BBTRANSLATION_H_*/
