#ifndef ARS_BBTRANSLATION_H_
#define ARS_BBTRANSLATION_H_

#include <ars/definitions.h>

#include <queue>

#include <Eigen/Dense>

namespace ars {

double distancePointBox(const Vector2& p,
                        const Vector2& boxMin,
                        const Vector2& boxMax);

struct Box {
    Vector2 min_;
    Vector2 max_;
    double lower_;
    double upper_;

    Box(const Vector2& min, const Vector2& max);

    Box(const Vector2& min,
        const Vector2& max,
        const VectorVector2& ptsSrc,
        const VectorVector2& ptsDst);

    virtual ~Box();

    void computeBounds(const VectorVector2& ptsSrc,
                       const VectorVector2& ptsDst);
};

std::ostream& operator<<(std::ostream& out, const Box& box);

class BBTranslation {
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
    void compute();

    void setTranslMinMax(const ars::Vector2& translMin,
                         const ars::Vector2& translMax);

    /**
     * @brief Set points src
     */
    void setPtsSrc(const ars::VectorVector2& pts);

    /**
     * @brief Set points dst
     */
    void setPtsDst(const ars::VectorVector2& pts);

    /**
     * @brief Set pts src and dst
     */
    void setPts(const ars::VectorVector2& ptsS, const ars::VectorVector2& ptsD);

   private:
    ars::Vector2 translMin_;
    ars::Vector2 translMax_;

    ars::VectorVector2 ptsSrc_;
    ars::VectorVector2 ptsDst_;
};
}  // namespace ars

#endif /*ARS_BBTRANSLATION_H_*/
