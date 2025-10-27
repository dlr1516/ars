#include <ars/FourierLowerUpperLut.h>
#include <ars/MortonSort.h>
#include <ars/definitions.h>
#include <ars/functions.h>

namespace ars {

FourierLowerUpperLut::FourierLowerUpperLut()
    : intervals_(), sinusoids_(), levelNum_(0) {}

FourierLowerUpperLut::FourierLowerUpperLut(const std::vector<double>& coeffs,
                                           size_t levelNum)
    : intervals_(), sinusoids_(), levelNum_(0) {
    init(coeffs, levelNum);
}

FourierLowerUpperLut::~FourierLowerUpperLut() {}

void FourierLowerUpperLut::init(const std::vector<double>& coeffs,
                                size_t levelNum) {
    auto comp = [](const CriticalPoint& cp1, const CriticalPoint& cp2) -> bool {
        return cp1.theta > cp2.theta;
    };
    std::priority_queue<CriticalPoint, std::vector<CriticalPoint>,
                        decltype(comp)>
        queue;
    Sinusoid ci;
    double period, yLowerInc, yUpperInc, yLowerDec, yUpperDec;

    // We expect that levelNum_ is not required in this implementation, but we
    // keep it for compatibility with the previous one.
    levelNum_ = levelNum;
    ARS_ASSERT_VAR2(coeffs.size() % 2 == 0, "Coeffs size must be even",
                    coeffs.size());
    size_t sinusoidNum = coeffs.size() / 2;

    // Fourier series: a0 + sum_i{ a_i*cos(2*i*x) + bi*sin(2*i*x) }
    // We convert it to: sum_i{ module_i*cos(2*i*x - phase_i) }
    ci.module = coeffs[0];
    ci.phase = 0.0;
    ci.phaseMin = 0.0;
    ci.order = 0;
    ci.increasing = true;
    sinusoids_.push_back(ci);
    for (size_t i = 1; i < sinusoidNum; ++i) {
        period = M_PI / i;

        ci.order = i;
        ci.module = std::sqrt(coeffs[2 * i] * coeffs[2 * i] +
                              coeffs[2 * i + 1] * coeffs[2 * i + 1]);
        // Phase normalized in [0, 2pi)
        ci.phase =
            std::fmod(std::atan2(coeffs[2 * i + 1], coeffs[2 * i]) + 2.0 * M_PI,
                      2.0 * M_PI);
        // PhaseMin is the first occurrence in interval [0, pi) of max of
        // sinsuoid with argument 2*i*x - phaseMin = 0 => x = phaseMin / (2*i)
        ci.phaseMin = std::fmod(ci.phase, M_PI / i);
        if (ci.phaseMin < 0.5 * period)
            ci.increasing = true;
        else
            ci.increasing = false;
        sinusoids_.push_back(ci);

        // Inserts the point where next change of monotonicity occurs
        // in the queue
        CriticalPoint cp;
        cp.order = i;
        cp.theta = ci.phaseMin;
        queue.push(cp);
    }

    // Process the queue to extract the critical points in [0, pi)
    double thetaPrev = 0.0;
    while (!queue.empty()) {
        CriticalPoint cp = queue.top();
        queue.pop();

        yLowerInc = 0.0;
        yUpperInc = 0.0;
        yLowerDec = 0.0;
        yUpperDec = 0.0;
        for (size_t k = 0; k < sinusoids_.size(); ++k) {
            if (sinusoids_[k].increasing) {
                yLowerInc += sinusoids_[k].module *
                             std::cos(2 * k * thetaPrev - sinusoids_[k].phase);
                yUpperInc += sinusoids_[k].module *
                             std::cos(2 * k * cp.theta - sinusoids_[k].phase);
            } else {
                yLowerDec += sinusoids_[k].module *
                             std::cos(2 * k * cp.theta - sinusoids_[k].phase);
                yUpperDec += sinusoids_[k].module *
                             std::cos(2 * k * thetaPrev - sinusoids_[k].phase);
            }

            if (k == cp.order) {
                // Update the sinusoid monotonicity
                sinusoids_[k].increasing = !sinusoids_[k].increasing;

                // Insert the next critical point in the queue
                CriticalPoint cpNext;
                cpNext.order = k;
                cpNext.theta = cp.theta + (0.5 * M_PI / k);
                if (cpNext.theta < M_PI) {
                    queue.push(cpNext);
                }
            }
        }

        Interval interval;
        interval.thetaMin = thetaPrev;
        interval.thetaMax = cp.theta;
        interval.yLower =
            std::min(yLowerInc + yLowerDec, yUpperInc + yUpperDec);
        interval.yUpper =
            std::max(yLowerInc + yLowerDec, yUpperInc + yUpperDec);
        intervals_.push_back(interval);
        ARS_VAR4(interval.thetaMin, interval.thetaMax, interval.yLower,
                 interval.yUpper);

        thetaPrev = cp.theta;
    }
}

void FourierLowerUpperLut::findLU(double xMin,
                                  double xMax,
                                  double& yLower,
                                  double& yUpper) const {}

// FourierLowerUpperLut::FourierLowerUpperLut()
//     : levelNum_(0), intervalNum_(0), dx_(0.0) {}

// FourierLowerUpperLut::FourierLowerUpperLut(const std::vector<double>& coeffs,
//                                            size_t levelNum) {
//     init(coeffs, levelNum);
// }

// FourierLowerUpperLut::~FourierLowerUpperLut() {}

// void FourierLowerUpperLut::init(const std::vector<double>& coeffs,
//                                 size_t levelNum) {
//     levelNum_ = levelNum;
//     intervalNum_ = 1 << levelNum_;
//     int treeSize = (1 << (levelNum_ + 1)) - 1;
//     dx_ = M_PI / static_cast<double>(intervalNum_);

//     luValues_.resize(intervalNum_);
//     intervals_.resize(treeSize);

//     size_t leafStart = levelStart(levelNum_);
//     for (int i = 0; i < intervalNum_; ++i) {
//         double xMin = i * dx_;
//         double xMax = (i + 1) * dx_;

//         double yLower, yUpper;
//         findLUFourier(coeffs, xMin, xMax, yLower, yUpper);

//         luValues_[i].lower = yLower;
//         luValues_[i].upper = yUpper;

//         intervals_[leafStart + i].idxL = i;
//         intervals_[leafStart + i].idxU = i;

//         size_t h = i;
//         size_t c = leafStart + i;
//         while (h > 0) {
//             size_t p = parent(c);
//             size_t cl = childLeft(p);
//             size_t cr = childRight(p);

//             size_t idxL = (luValues_[intervals_[cl].idxL].lower <
//                            luValues_[intervals_[cr].idxL].lower)
//                               ? intervals_[cl].idxL
//                               : intervals_[cr].idxL;

//             size_t idxU = (luValues_[intervals_[cl].idxU].upper >
//                            luValues_[intervals_[cr].idxU].upper)
//                               ? intervals_[cl].idxU
//                               : intervals_[cr].idxU;

//             intervals_[p].idxL = idxL;
//             intervals_[p].idxU = idxU;

//             h = h >> 1;
//             c = p;
//         }
//     }

//     // Debug print:
//     for (size_t i = 0; i < luValues_.size(); ++i) {
//         size_t idxLower = intervals_[i].idxL;
//         size_t idxUpper = intervals_[i].idxU;
//         double lower = luValues_[idxLower].lower;
//         double upper = luValues_[idxUpper].upper;
//         ARS_PRINT("node " << i << " idxLower " << idxLower << ", idxUpper "
//                           << idxUpper << "[" << lower << "," << upper <<
//                           "]");
//     }
// }

// void FourierLowerUpperLut::findLU(double xMin,
//                                   double xMax,
//                                   double& yLower,
//                                   double& yUpper) const {
//     if (xMin > xMax)
//         std::swap(xMin, xMax);

//     if (xMax - xMin >= M_PI) {
//         size_t idxL = intervals_[0].idxL;
//         size_t idxU = intervals_[0].idxU;
//         yLower = luValues_[idxL].lower;
//         yUpper = luValues_[idxU].upper;
//         return;
//     } else {
//         size_t idxMin = (intervalNum_ + (int)floor(xMin / dx_)) %
//         intervalNum_; size_t idxMax = (intervalNum_ + (int)ceil(xMax / dx_))
//         % intervalNum_;

//         ARS_VAR4(xMin, xMax, idxMin, idxMax);
//         std::cout << "xMin " << xMin << " in [" << (idxMin * dx_) << ", "
//                   << (idxMin + 1) * dx_ << "], xMax " << xMax << " in ["
//                   << (idxMax - 1) * dx_ << ", " << idxMax * dx_ << "]"
//                   << std::endl;
//         if (idxMin <= idxMax) {
//             findLUTree(idxMin, idxMax, yLower, yUpper);
//         } else {
//             size_t ancestorL = findCommonAncestor(idxMin, intervalNum_ - 1);
//             size_t ancestorU = findCommonAncestor(0, idxMax - 1);

//             double yLower1, yUpper1, yLower2, yUpper2;
//             findLUTree(idxMin, intervalNum_ - 1, yLower1, yUpper1);
//             findLUTree(0, idxMax, yLower2, yUpper2);

//             yLower = std::min(yLower1, yLower2);
//             yUpper = std::max(yUpper1, yUpper2);
//         }
//     }
// }

// void FourierLowerUpperLut::findLUTree(size_t idxMin,
//                                       size_t idxMax,
//                                       double& lower,
//                                       double& upper) const {
//     size_t idxLow, idxMid, idxUpp;
//     size_t ancestor, idxL, idxU;
//     double lower1, upper1, lower2, upper2;

//     intervalPow2(idxMin, idxMax, idxLow, idxMid, idxUpp);

//     ARS_VAR5(idxMin, idxMax, idxLow, idxMid, idxUpp);
//     std::cout << " idxMin:  " << std::bitset<32>(idxMin) << "\n"
//               << " idxMax:  " << std::bitset<32>(idxMax) << "\n"
//               << " idxLow:  " << std::bitset<32>(idxLow) << "\n"
//               << " idxMid:  " << std::bitset<32>(idxMid) << "\n"
//               << " idxUpp:  " << std::bitset<32>(idxUpp) << "\n"
//               << std::endl;

//     if (idxLow == idxMin && idxUpp == idxMax) {
//         ancestor = findCommonAncestor(idxMin, idxMax - 1);
//         idxL = intervals_[ancestor].idxL;
//         idxU = intervals_[ancestor].idxU;
//         lower = luValues_[idxL].lower;
//         upper = luValues_[idxU].upper;
//         return;
//     }

//     findLUTree(idxMin, idxMid - 1, lower1, upper1);
//     findLUTree(idxMid, idxMax, lower2, upper2);

//     lower = std::min(lower1, lower2);
//     upper = std::max(upper1, upper2);
// }

// size_t FourierLowerUpperLut::findCommonAncestor(size_t idxL,
//                                                 size_t idxU) const {
//     size_t nodeL = levelStart(levelNum_) + idxL;
//     size_t nodeU = levelStart(levelNum_) + idxU;

//     while (nodeL != nodeU) {
//         if (nodeL > nodeU) {
//             nodeL = parent(nodeL);
//         } else {
//             nodeU = parent(nodeU);
//         }
//     }

//     return nodeL;
// }

}  // namespace ars