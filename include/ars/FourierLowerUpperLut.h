#ifndef ARS_LU_BOUNDS_LUT_H_
#define ARS_LU_BOUNDS_LUT_H_

#include <ars/definitions.h>
#include <functional>
#include <iostream>
#include <vector>

namespace ars {

class FourierLowerUpperLut {
   public:
    /**
     * Struct Sinusoid stores the parameters of a sinusoid term in the Fourier
     * series.
     */
    struct Sinusoid {
        size_t order;
        double module;
        double phase;
        bool increasing;
    };
    using Sinusoids = std::vector<Sinusoid>;

    /**
     * Struct PointInterval stores the point where the monotonicity of
     * sinusoid of order k changes.
     */
    struct CriticalPoint {
        size_t order;
        double x;
    };

    /**
     * Struct Interval stores the lower and upper bounds of the Fourier series
     * in an interval where all the sinusoid terms are monotonic.
     */
    struct Interval {
        double xMin;
        double xMax;
        double yLower;
        double yUpper;
    };
    using Intervals = std::vector<Interval>;

    struct IndexNode {
        double xMin;
        double xMax;
        size_t idxYL;
        size_t idxYU;
        struct IndexNode* left;
        struct IndexNode* right;
    };

    /**
     * Default constructor.
     */
    FourierLowerUpperLut();

    /**
     * Constructor with Fourier coefficients.
     * @param coeffs Fourier coefficients
     */
    FourierLowerUpperLut(const std::vector<double>& coeffs);

    /**
     * @brief Destroy the Fourier Lower Upper Lut object
     */
    ~FourierLowerUpperLut();

    /**
     * @brief Initialize the LUT with Fourier coefficients
     */
    void init(const std::vector<double>& coeffs);

    /**
     * @brief Find the lower and upper bounds for the given x interval
     *
     * @param xMin lower bound of x
     * @param xMax upper bound of x
     * @param yLower output lower bound of y
     * @param yUpper output upper bound of y
     */
    void findLU(double xMin, double xMax, double& yLower, double& yUpper) const;

    /**
     * @brief Added for debug only. It will be removed!
     *
     * @return const std::vector<Interval>&
     */
    const std::vector<Interval>& intervals() const { return intervals_; }

    /**
     * @brief Export the LUT data for plotting
     *
     * @param out output stream
     */
    void exportPlot(std::ostream& out);

   protected:
    Intervals intervals_;
    Sinusoids sinusoids_;
    IndexNode* tree_;

    IndexNode* buildTree(size_t idxBeg, size_t idxEnd);

    void removeTree(IndexNode* node);

    void findLUTree(IndexNode* node,
                    double xMin,
                    double xMax,
                    double& yLower,
                    double& yUpper) const;
};

}  // namespace ars

#endif