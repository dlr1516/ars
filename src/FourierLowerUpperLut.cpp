#include <ars/FourierLowerUpperLut.h>
#include <ars/MortonSort.h>
#include <ars/functions.h>

namespace ars {

FourierLowerUpperLut::FourierLowerUpperLut()
    : levelNum_(0), intervalNum_(0), dx_(0.0) {}

FourierLowerUpperLut::FourierLowerUpperLut(const std::vector<double>& coeffs,
                                           size_t levelNum) {
    init(coeffs, levelNum);
}

FourierLowerUpperLut::~FourierLowerUpperLut() {}

void FourierLowerUpperLut::init(const std::vector<double>& coeffs,
                                size_t levelNum) {
    levelNum_ = levelNum;
    intervalNum_ = 1 << levelNum_;
    int treeSize = (1 << (levelNum_ + 1)) - 1;
    dx_ = M_PI / static_cast<double>(intervalNum_);

    luValues_.resize(intervalNum_);
    intervals_.resize(treeSize);

    size_t leafStart = levelStart(levelNum_);
    for (int i = 0; i < intervalNum_; ++i) {
        double xMin = i * dx_;
        double xMax = (i + 1) * dx_;

        double yLower, yUpper;
        findLUFourier(coeffs, xMin, xMax, yLower, yUpper);

        luValues_[i].lower = yLower;
        luValues_[i].upper = yUpper;

        intervals_[leafStart + i].idxL = i;
        intervals_[leafStart + i].idxU = i;

        size_t h = i;
        size_t c = leafStart + i;
        while (h > 0) {
            size_t p = parent(c);
            size_t cl = childLeft(p);
            size_t cr = childRight(p);

            size_t idxL = (luValues_[intervals_[cl].idxL].lower <
                           luValues_[intervals_[cr].idxL].lower)
                              ? intervals_[cl].idxL
                              : intervals_[cr].idxL;

            size_t idxU = (luValues_[intervals_[cl].idxU].upper >
                           luValues_[intervals_[cr].idxU].upper)
                              ? intervals_[cl].idxU
                              : intervals_[cr].idxU;

            intervals_[p].idxL = idxL;
            intervals_[p].idxU = idxU;

            h = h >> 1;
            c = p;
        }
    }

    // Debug print:
    for (size_t i = 0; i < luValues_.size(); ++i) {
        size_t idxLower = intervals_[i].idxL;
        size_t idxUpper = intervals_[i].idxU;
        double lower = luValues_[idxLower].lower;
        double upper = luValues_[idxUpper].upper;
        ARS_PRINT("node " << i << " idxLower " << idxLower << ", idxUpper "
                          << idxUpper << "[" << lower << "," << upper << "]");
    }
}

void FourierLowerUpperLut::findLU(double xMin,
                                  double xMax,
                                  double& yLower,
                                  double& yUpper) const {
    if (xMin > xMax)
        std::swap(xMin, xMax);

    if (xMax - xMin >= M_PI) {
        size_t idxL = intervals_[0].idxL;
        size_t idxU = intervals_[0].idxU;
        yLower = luValues_[idxL].lower;
        yUpper = luValues_[idxU].upper;
        return;
    } else {
        size_t idxMin = (intervalNum_ + (int)floor(xMin / dx_)) % intervalNum_;
        size_t idxMax = (intervalNum_ + (int)ceil(xMax / dx_)) % intervalNum_;

        if (idxMin <= idxMax) {
            findLUTree(idxMin, idxMax, yLower, yUpper);
        } else {
            size_t ancestorL = findCommonAncestor(idxMin, intervalNum_ - 1);
            size_t ancestorU = findCommonAncestor(0, idxMax - 1);

            double yLower1, yUpper1, yLower2, yUpper2;
            findLUTree(idxMin, intervalNum_ - 1, yLower1, yUpper1);
            findLUTree(0, idxMax, yLower2, yUpper2);

            yLower = std::min(yLower1, yLower2);
            yUpper = std::max(yUpper1, yUpper2);
        }
    }
}

void FourierLowerUpperLut::findLUTree(size_t idxMin,
                                      size_t idxMax,
                                      double& lower,
                                      double& upper) const {
    size_t idxLow, idxMid, idxUpp;
    size_t ancestor, idxL, idxU;
    double lower1, upper1, lower2, upper2;

    intervalPow2(idxMin, idxMax, idxLow, idxMid, idxUpp);

    if (idxLow == idxMin && idxUpp == idxMax) {
        ancestor = findCommonAncestor(idxMin, idxMax - 1);
        idxL = intervals_[ancestor].idxL;
        idxU = intervals_[ancestor].idxU;
        lower = luValues_[idxL].lower;
        upper = luValues_[idxU].upper;
        return;
    }

    findLU(idxMin, idxMid, lower1, upper1);
    findLU(idxMid, idxMax, lower2, upper2);

    lower = std::min(lower1, lower2);
    upper = std::max(upper1, upper2);
}

size_t FourierLowerUpperLut::findCommonAncestor(size_t idxL,
                                                size_t idxU) const {
    size_t nodeL = levelStart(levelNum_) + idxL;
    size_t nodeU = levelStart(levelNum_) + idxU;

    while (nodeL != nodeU) {
        if (nodeL > nodeU) {
            nodeL = parent(nodeL);
        } else {
            nodeU = parent(nodeU);
        }
    }

    return nodeL;
}

}  // namespace ars