#ifndef ARS_BBTRANSLATION_H_
#define ARS_BBTRANSLATION_H_

#include <ars/definitions.h>

#include <queue>

#include <Eigen/Dense>
#include <Eigen/Core>

namespace ars
{

    double distancePointBox(const Vector2 &p,
                            const Vector2 &boxMin,
                            const Vector2 &boxMax);

    struct Box
    {
        Vector2 min_;
        Vector2 max_;
        double lower_;
        double upper_;
        double eps_;

        Box(const Vector2 &min, const Vector2 &max, const double eps);

        Box(const Vector2 &min,
            const Vector2 &max,
            const VectorVector2 &ptsSrc,
            const VectorVector2 &ptsDst,
            const double eps);

        virtual ~Box();

        void computeBoundsNaive(const VectorVector2 &ptsSrc,
                                const VectorVector2 &ptsDst);

        void computeBoundsInlier(const VectorVector2 &ptsSrc,
                                 const VectorVector2 &ptsDst);
    };

    std::ostream &operator<<(std::ostream &out, const Box &box);

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
