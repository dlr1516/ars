#ifndef HISTOGRAM_CIRCULAR_CORRELATION_H
#define HISTOGRAM_CIRCULAR_CORRELATION_H

#include <iostream>
#include <Eigen/Dense>

namespace ars {

    class HistogramCircularCorrelation {
    public:
        /** Constructor. 
         */
        HistogramCircularCorrelation();

        /** Computes the value of histogram correlation.
         */
        void computeHistogramCorrelation(const Eigen::VectorXd& h1, const Eigen::VectorXd& h2, int shift, double& corr);

        /** Computes the histogram shift that maximize histogram correlation.
         */
        void computeShiftBruteForce(const Eigen::VectorXd& h1, const Eigen::VectorXd& h2, int shiftMax, int& shift, double& corr);

        /** Computes the value of histogram correlation.
         */
        template <typename Hist, typename Value>
        void computeHistSimpleCorr(const Hist& h1, const Hist& h2, int shift, Value& corr) {
            assert(h1.size() == h2.size());
            int n = h1.size();
            int beg = std::max(0, shift);
            int end = std::min(n, n + shift);
            corr = 0;
            //    std::cout << __FILE__ << ": shift " << shift << " n " << n 
            //      << ": corr = h1(" << beg << ")*h2(" << (beg-shift) << ") "
            //      << "+ ... + h1(" << (end-1) << ")*h2(" << (end-1-shift) << ") " << std::endl;
            for (int i = beg; i < end; ++i) {
                corr += h1(i) * h2(i - shift);
            }
        }

        /** Computes the histogram shift that maximize histogram correlation.
         */
        template <typename Hist, typename Value>
        void computeHistSimpleShift(const Hist& h1, const Hist& h2, int shiftMax, int& shift, Value& corr) {
            assert(h1.size() == h2.size());
            Value corrTmp = 0;
            corr = 0;
            std::cout << __FILE__ << "," << __LINE__ << ": find histogram correlation" << std::endl;
            for (int i = -shiftMax; i < shiftMax; ++i) {
                computeHistSimpleCorr(h1, h2, i, corrTmp);
                if (corrTmp > corr) {
                    corr = corrTmp;
                    shift = i;
                }
                std::cout << "  shift " << i << " corr " << corrTmp << " best shift " << shift << " best corr " << corr << std::endl;
            }
        }

        /** Computes the value of histogram correlation.
         */
        template <typename Hist, typename Value>
        void computeHistCircularCorr(const Hist& h1, const Hist& h2, int shift, Value& corr) {
            assert(h1.size() == h2.size());
            int n = h1.size();
            shift = shift % n;
            corr = 0;
            for (int i = 0; i < n; ++i) {
                corr += h1(i) * h2((i + shift + n) % n);
            }
        }

        /** Computes the histogram shift that maximize histogram correlation.
         */
        template <typename Hist, typename Value>
        void computeHistCircularShift(const Hist& h1, const Hist& h2, int shiftMax, int& shift, Value& corr) {
            assert(h1.size() == h2.size());
            Value corrTmp = 0;
            corr = 0;
            for (int i = -shiftMax; i < shiftMax; ++i) {
                computeHistCircularCorr(h1, h2, i, corrTmp);
                if (corrTmp > corr) {
                    corr = corrTmp;
                    shift = i;
                }
            }
        }
    };

} // end of namespace

#endif
