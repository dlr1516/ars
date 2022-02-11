#include <ars/HistogramCircularCorrelation.h>

namespace cuars {

    HistogramCircularCorrelation::HistogramCircularCorrelation() {
    }

    void HistogramCircularCorrelation::computeHistogramCorrelation(const Eigen::VectorXd& h1, const Eigen::VectorXd& h2, int shift, double& corr) {
        assert(h1.size() == h2.size());
        int n = h1.size();
        shift = shift % n;
        corr = 0.0;
        for (int i = 0; i < n; ++i) {
            corr += h1(i) * h2((i + shift + n) % n);
        }
    }

    void HistogramCircularCorrelation::computeShiftBruteForce(const Eigen::VectorXd& h1, const Eigen::VectorXd& h2, int shiftMax, int& shift, double& corr) {
        assert(h1.size() == h2.size());
        double corrTmp = 0.0;
        corr = 0.0;
        shift = 0;
        for (int i = -shiftMax; i < shiftMax; ++i) {
            computeHistogramCorrelation(h1, h2, i, corrTmp);
            if (corrTmp > corr) {
                corr = corrTmp;
                shift = i;
            }
        }
    }

} // end of namespace
