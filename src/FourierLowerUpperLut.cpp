#include <ars/FourierLowerUpperLut.h>
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
}

}  // namespace ars