#ifndef ARS_LU_BOUNDS_LUT_H_
#define ARS_LU_BOUNDS_LUT_H_

#include <functional>
#include <iostream>
#include <vector>

namespace ars {

class FourierLowerUpperLut {
   public:
    struct LUIndex {
        size_t idxL;
        size_t idxU;
    };

    struct LUValues {
        double lower;
        double upper;
    };

    FourierLowerUpperLut();

    FourierLowerUpperLut(const std::vector<double>& coeffs, size_t levelNum);

    ~FourierLowerUpperLut();

    void init(const std::vector<double>& coeffs, size_t levelNum);

    void findLU(double xMin, double xMax, double& yLower, double& yUpper) const;

   private:
    std::vector<LUValues> luValues_;
    std::vector<LUIndex> intervals_;
    size_t levelNum_;
    size_t intervalNum_;
    double dx_;

    void findLUTree(size_t idxMin,
                    size_t idxMax,
                    double& lower,
                    double& upper) const;

    /**
     * The vector intervals_ stores a tree organized into levels.
     * The left and right children of a node p are:
     *    childLeft(p) = 2*p + 1
     *    childRight(p) = 2*p + 2
     * Thus the parent of a child c is:
     *   parent(c) = (c - 1) / 2
     */
    inline size_t parent(size_t c) const { return ((c - 1) >> 1); }

    inline size_t childLeft(size_t p) const { return ((p << 1) + 1); }

    inline size_t childRight(size_t p) const { return ((p << 1) + 2); }

    inline size_t levelStart(size_t level) const { return (1 << level) - 1; }

    size_t findCommonAncestor(size_t idxL, size_t idxU) const;
};

}  // namespace ars

#endif