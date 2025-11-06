#include <ars/FourierLowerUpperLut.h>
#include <ars/MortonSort.h>
#include <ars/functions.h>

#define RAD2DEG(X) (180.0 / M_PI * (X))

namespace ars {

FourierLowerUpperLut::FourierLowerUpperLut()
    : intervals_(), sinusoids_(), tree_(nullptr) {}

FourierLowerUpperLut::FourierLowerUpperLut(const std::vector<double>& coeffs)
    : intervals_(), sinusoids_(), tree_(nullptr) {
    init(coeffs);
}

FourierLowerUpperLut::~FourierLowerUpperLut() {
    removeTree(tree_);
}

void FourierLowerUpperLut::init(const std::vector<double>& coeffs) {
    auto comp = [](const CriticalPoint& cp1, const CriticalPoint& cp2) -> bool {
        return cp1.x > cp2.x;
    };
    std::priority_queue<CriticalPoint, std::vector<CriticalPoint>,
                        decltype(comp)>
        queue;
    Sinusoid ci;
    double period, xBeg, xPrev;
    double yLowerInc, yUpperInc, yLowerDec, yUpperDec, y1, y2;

    // We expect that levelNum_ is not required in this implementation, but we
    // keep it for compatibility with the previous one.
    ARS_ASSERT_VAR2(coeffs.size() % 2 == 0, "Coeffs size must be even",
                    coeffs.size());
    size_t sinusoidNum = coeffs.size() / 2;

    // Fourier series: a0 + sum_i{ a_i*cos(2*i*x) + bi*sin(2*i*x) }
    // We convert it to: sum_i{ module_i*cos(2*i*x - phase_i) }
    ci.module = coeffs[0];
    ci.phase = 0.0;
    // ci.phaseMin = 0.0;
    ci.order = 0;
    ci.increasing = true;
    sinusoids_.push_back(ci);
    for (size_t k = 1; k < sinusoidNum; ++k) {
        period = M_PI / k;

        ci.order = k;
        ci.module = std::sqrt(coeffs[2 * k] * coeffs[2 * k] +
                              coeffs[2 * k + 1] * coeffs[2 * k + 1]);
        // Phase normalized in [0, 2pi)
        ci.phase =
            std::fmod(std::atan2(coeffs[2 * k + 1], coeffs[2 * k]) + 2.0 * M_PI,
                      2.0 * M_PI);
        // xBeg is the first occurrence in interval [0, pi) of min or max
        // of sinsuoid with argument 2*k*x - phase
        //   2 * k * xBeg - phase = j * pi   (j is 0 or 1)
        xBeg = ci.phase / (2.0 * k);
        if (xBeg <= 0.5 * period)
            ci.increasing = true;
        else {
            ci.increasing = false;
            xBeg -= 0.5 * period;
        }

        sinusoids_.push_back(ci);
        ARS_VAR5(k, ci.module, RAD2DEG(ci.phase), RAD2DEG(xBeg),
                 RAD2DEG(period));

        // Inserts the point where next change of monotonicity occurs
        // in the queue
        CriticalPoint cp;
        cp.order = k;
        cp.x = xBeg;
        queue.push(cp);
    }

    // Process the queue to extract the critical points in [0, pi)
    xPrev = 0.0;
    while (!queue.empty()) {
        CriticalPoint cp = queue.top();
        queue.pop();

        yLowerInc = 0.0;
        yUpperInc = 0.0;
        yLowerDec = 0.0;
        yUpperDec = 0.0;
        for (size_t k = 0; k < sinusoids_.size(); ++k) {
            y1 = sinusoids_[k].module *
                 std::cos(2 * k * xPrev - sinusoids_[k].phase);
            y2 = sinusoids_[k].module *
                 std::cos(2 * k * cp.x - sinusoids_[k].phase);

            ARS_VAR6(k, sinusoids_[k].increasing, RAD2DEG(xPrev), RAD2DEG(cp.x),
                     y1, y2);

            if (sinusoids_[k].increasing) {
                yLowerInc += y1;
                yUpperInc += y2;
            } else {
                yLowerDec += y2;
                yUpperDec += y1;
            }

            if (k == cp.order) {
                // Update the sinusoid monotonicity
                sinusoids_[k].increasing = !sinusoids_[k].increasing;

                // Insert the next critical point in the queue
                CriticalPoint cpNext;
                cpNext.order = k;
                cpNext.x = cp.x + (0.5 * M_PI / k);
                if (cpNext.x < M_PI) {
                    queue.push(cpNext);
                }
            }
        }

        Interval interval;
        interval.xMin = xPrev;
        interval.xMax = cp.x;
        interval.yLower =
            std::min(yLowerInc + yLowerDec, yUpperInc + yUpperDec);
        interval.yUpper =
            std::max(yLowerInc + yLowerDec, yUpperInc + yUpperDec);
        intervals_.push_back(interval);
        ARS_VAR4(interval.xMin, interval.xMax, interval.yLower,
                 interval.yUpper);

        xPrev = cp.x;
    }

    // Final interval from last critical point to pi
    yLowerInc = 0.0;
    yUpperInc = 0.0;
    yLowerDec = 0.0;
    yUpperDec = 0.0;
    for (size_t k = 0; k < sinusoids_.size(); ++k) {
        y1 = sinusoids_[k].module *
             std::cos(2 * k * xPrev - sinusoids_[k].phase);
        y2 =
            sinusoids_[k].module * std::cos(2 * k * M_PI - sinusoids_[k].phase);

        ARS_VAR6(k, sinusoids_[k].increasing, RAD2DEG(xPrev), RAD2DEG(M_PI), y1,
                 y2);

        if (sinusoids_[k].increasing) {
            yLowerInc += y1;
            yUpperInc += y2;
        } else {
            yLowerDec += y2;
            yUpperDec += y1;
        }
    }
    Interval interval;
    interval.xMin = xPrev;
    interval.xMax = M_PI;
    interval.yLower = std::min(yLowerInc + yLowerDec, yUpperInc + yUpperDec);
    interval.yUpper = std::max(yLowerInc + yLowerDec, yUpperInc + yUpperDec);
    intervals_.push_back(interval);
    ARS_VAR4(interval.xMin, interval.xMax, interval.yLower, interval.yUpper);

    if (tree_ != nullptr) {
        removeTree(tree_);
    }
    tree_ = buildTree(0, intervals_.size());
}

void FourierLowerUpperLut::findLU(double xMin,
                                  double xMax,
                                  double& yLower,
                                  double& yUpper) const {
    ARS_ASSERT(tree_ != nullptr);
    yLower = intervals_[tree_->idxYL].yLower;
    yUpper = intervals_[tree_->idxYU].yUpper;
    findLUTree(tree_, xMin, xMax, yLower, yUpper);
}

void FourierLowerUpperLut::exportPlot(std::ostream& out) {
    const size_t NUM = 720;
    double dx = M_PI / static_cast<double>(NUM);

    for (size_t k = 0; k < sinusoids_.size(); ++k) {
        out << "set term wxt " << k << "\n";
        out << "set title 'Sinusoid order " << sinusoids_[k].order << "'\n";
        out << "plot '-' title 'sinusoid' w l\n";
        for (size_t i = 0; i <= NUM; ++i) {
            double x = i * dx;
            double y =
                sinusoids_[k].module *
                std::cos(2 * sinusoids_[k].order * x - sinusoids_[k].phase);
            out << RAD2DEG(x) << " " << y << "\n";
        }
        out << "e\n";
    }
}

void FourierLowerUpperLut::printTree(std::ostream& out) const {
    ARS_PRINT("Interval Index Tree:");
    printTree(out, tree_, 0);
}

// ------------------------------------------------------------------
// PRIVATE METHODS
// ------------------------------------------------------------------

FourierLowerUpperLut::IndexNode* FourierLowerUpperLut::buildTree(
    size_t idxBeg,
    size_t idxEnd) {
    IndexNode* node;
    size_t idxMid;

    if (idxBeg == idxEnd) {
        return nullptr;
    } else if (idxEnd == idxBeg + 1) {
        node = new IndexNode;
        node->xMin = intervals_[idxBeg].xMin;
        node->xMax = intervals_[idxBeg].xMax;
        node->idxYL = idxBeg;
        node->idxYU = idxBeg;
        node->left = nullptr;
        node->right = nullptr;
        return node;
    } else {
        node = new IndexNode;
        idxMid = (idxBeg + idxEnd) / 2 + (idxBeg + idxEnd) % 2;
        node->left = buildTree(idxBeg, idxMid);
        node->right = buildTree(idxMid, idxEnd);

        ARS_ASSERT(node->left != nullptr || node->right != nullptr);

        node->xMin = node->left->xMin;
        node->xMax = node->right->xMax;
        node->idxYL = intervals_[node->left->idxYL].yLower <
                              intervals_[node->right->idxYL].yLower
                          ? node->left->idxYL
                          : node->right->idxYL;
        node->idxYU = intervals_[node->left->idxYU].yUpper >
                              intervals_[node->right->idxYU].yUpper
                          ? node->left->idxYU
                          : node->right->idxYU;
        return node;
    }
    return nullptr;
}

void FourierLowerUpperLut::removeTree(IndexNode* node) {
    if (node != nullptr) {
        removeTree(node->left);
        removeTree(node->right);
        delete node;
    }
}

void FourierLowerUpperLut::printTree(std::ostream& out,
                                     IndexNode* node,
                                     int level) const {
    if (node != nullptr) {
        out << std::string(level * 2, ' ') << "Node [" << RAD2DEG(node->xMin)
            << ", " << RAD2DEG(node->xMax) << "]: idxYL=" << node->idxYL << " ("
            << intervals_[node->idxYL].yLower << ")"
            << ", idxYU=" << node->idxYU << " ("
            << intervals_[node->idxYU].yUpper << ")" << std::endl;
        printTree(out, node->left, level + 1);
        printTree(out, node->right, level + 1);
    }
}

void FourierLowerUpperLut::findLUTree(IndexNode* node,
                                      double xMin,
                                      double xMax,
                                      double& yLower,
                                      double& yUpper) const {
    if (node == nullptr) {
        return;
    }
    if (xMax < node->xMin || xMin > node->xMax) {
        // No overlap
        return;
    } else if (xMin <= node->xMin && node->xMax <= xMax) {
        // Full overlap
        yLower = std::min(yLower, intervals_[node->idxYL].yLower);
        yUpper = std::max(yUpper, intervals_[node->idxYU].yUpper);
        ARS_PRINT("Full overlap node [" << RAD2DEG(node->xMin) << ", "
                                        << RAD2DEG(node->xMax)
                                        << "]: " << yLower << ", " << yUpper);
    } else {
        // Partial overlap
        findLUTree(node->left, xMin, xMax, yLower, yUpper);
        findLUTree(node->right, xMin, xMax, yLower, yUpper);
    }
}

}  // namespace ars