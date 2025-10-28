#ifndef ARS_LU_BOUNDS_LUT_H_
#define ARS_LU_BOUNDS_LUT_H_

#include <functional>
#include <iostream>
#include <vector>

namespace ars {

class FourierLowerUpperLut {
   public:
    struct Sinusoid {
        size_t order;
        double module;
        double phase;
        // double phaseMin;
        bool increasing;
    };

    /**
     * PointInterval stores the point where the monotonicity of
     * sinusoid of order k changes.
     */
    struct CriticalPoint {
        size_t order;
        double theta;
    };

    /**
     *
     */
    struct Interval {
        double thetaMin;
        double thetaMax;
        double yLower;
        double yUpper;
    };
    using Intervals = std::vector<Interval>;

    FourierLowerUpperLut();

    FourierLowerUpperLut(const std::vector<double>& coeffs, size_t levelNum);

    ~FourierLowerUpperLut();

    void init(const std::vector<double>& coeffs, size_t levelNum);

    void findLU(double xMin, double xMax, double& yLower, double& yUpper) const;

    /**
     * @brief Added for debug only. It will be removed!
     *
     * @return const std::vector<Interval>&
     */
    const std::vector<Interval>& intervals() const { return intervals_; }

    void exportPlot(std::ostream& out);

   private:
    Intervals intervals_;
    std::vector<Sinusoid> sinusoids_;
    size_t levelNum_;
};

// class FourierLowerUpperLut {
//    public:
//     struct LUIndex {
//         size_t idxL;
//         size_t idxU;
//     };

//     struct LUValues {
//         double lower;
//         double upper;
//     };

//     FourierLowerUpperLut();

//     FourierLowerUpperLut(const std::vector<double>& coeffs, size_t levelNum);

//     ~FourierLowerUpperLut();

//     void init(const std::vector<double>& coeffs, size_t levelNum);

//     void findLU(double xMin, double xMax, double& yLower, double& yUpper)
//     const;

//    private:
//     std::vector<LUValues> luValues_;
//     std::vector<LUIndex> intervals_;
//     size_t levelNum_;
//     size_t intervalNum_;
//     double dx_;

//     void findLUTree(size_t idxMin,
//                     size_t idxMax,
//                     double& lower,
//                     double& upper) const;

//     /**
//      * The vector intervals_ stores a tree organized into levels.
//      * The left and right children of a node p are:
//      *    childLeft(p) = 2*p + 1
//      *    childRight(p) = 2*p + 2
//      * Thus the parent of a child c is:
//      *   parent(c) = (c - 1) / 2
//      */
//     inline size_t parent(size_t c) const { return ((c - 1) >> 1); }

//     inline size_t childLeft(size_t p) const { return ((p << 1) + 1); }

//     inline size_t childRight(size_t p) const { return ((p << 1) + 2); }

//     inline size_t levelStart(size_t level) const { return (1 << level) - 1; }

//     size_t findCommonAncestor(size_t idxL, size_t idxU) const;
// };

}  // namespace ars

#endif