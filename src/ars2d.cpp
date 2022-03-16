/**
 * ARS - Angular Radon Spectrum 
 * Copyright (C) 2017 Dario Lodi Rizzini.
 *
 * ARS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * ARS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with ARS.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <ars/ars2d.h>
#include <boost/numeric/interval.hpp>

#include <ars/Profiler.h>


namespace cuars {

    // --------------------------------------------------------
    // ARSF FUNCTIONS
    // --------------------------------------------------------

    /** Computes the contribution of term (lambda,phi) to the Fourier coefficients of ARS
     * computed with downward recursion of modified Bessel functions.
     */
    void updateARSF2CoeffRecursDown(double lambda, double cth2, double sth2, double factor, int n, std::vector<double>& coeffs) {
        std::vector<double> pnebis(n + 1);
        double sgn, cth, sth, ctmp, stmp;

        // Fourier Coefficients 
        if (coeffs.size() != 2 * n + 2) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid size of Fourier coefficients vector " << coeffs.size() << " should be " << (2 * n + 2) << std::endl;
            return;
        }

        evaluatePnebiVector(n, lambda, pnebis);
        //ARS_PRINT(pnebis[0]);
        coeffs[0] += 0.5 * factor * pnebis[0];
        std::cout << "coeff0" << 0.5 * factor * pnebis[0] << std::endl;
        sgn = -1.0;
        cth = cth2;
        sth = sth2;
        for (int k = 1; k <= n; ++k) {
            std::cout << "coeff" << 2 * k << " " << factor * pnebis[k] * sgn * cth << std::endl;
            std::cout << "coeff" << 2 * k + 1 << " " << factor * pnebis[k] * sgn * cth << std::endl;
            coeffs[2 * k] += factor * pnebis[k] * sgn * cth;
            coeffs[2 * k + 1] += factor * pnebis[k] * sgn * sth;
            sgn = -sgn;
            ctmp = cth2 * cth - sth2 * sth;
            stmp = sth2 * cth + cth2 * sth;
            cth = ctmp;
            sth = stmp;
        }
    }

    void updateARSF2CoeffRecursDownLUT(double lambda, double cth2, double sth2, double factor, int n, const PnebiLUT& pnebiLUT, std::vector<double>& coeffs) {
        std::vector<double> pnebis(n + 1);
        double sgn, cth, sth, ctmp, stmp;

        // Fourier Coefficients 
        if (coeffs.size() != 2 * n + 2 || pnebiLUT.getOrderMax() < n) {
            std::cerr << __FILE__ << "," << __LINE__ << ": one of these conditions failed:"
                    << "\n  size of Fourier coefficients vector " << coeffs.size() << " should be " << (2 * n + 2)
                    << "\n  LUT max order is " << pnebiLUT.getOrderMax() << " >= " << n
                    << std::endl;
            return;
        }

        pnebiLUT.eval(lambda, pnebis);
        //ARS_PRINT(pnebis[0]);
        coeffs[0] += 0.5 * factor * pnebis[0];
        sgn = -1.0;
        cth = cth2;
        sth = sth2;
        for (int k = 1; k <= n; ++k) {
            coeffs[2 * k] += pnebis[k] * factor * sgn * cth;
            coeffs[2 * k + 1] += pnebis[k] * factor * sgn * sth;
            sgn = -sgn;
            ctmp = cth2 * cth - sth2 * sth;
            stmp = sth2 * cth + cth2 * sth;
            cth = ctmp;
            sth = stmp;
        }
    }

    void computeFourierCorr(const std::vector<double>& fourierSrc, const std::vector<double>& fourierDst, std::vector<double>& fourierCor) {
        int n;

        if (fourierSrc.size() % 2 != 0 || fourierDst.size() % 2 != 0 || fourierSrc.size() != fourierDst.size()) {
            std::cerr << __FILE__ << "," << __LINE__ << ": the number of ARSF coefficients must be even and equal: fourierSrc " << fourierSrc.size()
                    << ", fourierDst " << fourierDst.size() << std::endl;
            return;
        }
        n = (fourierSrc.size() / 2) - 1;

        // Resizes the correlation coefficients, if required
        if (fourierCor.size() != 2 * n + 2) {
            fourierCor.resize(2 * n + 2);
        }

        // Computes the coefficients
        for (int k = 0; k <= n; ++k) {
            fourierCor[2 * k] = 0.5 * (fourierSrc[2 * k] * fourierDst[2 * k] + fourierSrc[2 * k + 1] * fourierDst[2 * k + 1]);
            fourierCor[2 * k + 1] = 0.5 * (fourierSrc[2 * k] * fourierDst[2 * k + 1] - fourierSrc[2 * k + 1] * fourierDst[2 * k]);
        }
    }

    // --------------------------------------------------------
    // ARS 2D CLASS
    // --------------------------------------------------------

    //    std::array<std::string, 2> const AngularRadonSpectrum2d::MODE_NAME{"PNEBI_DOWNWARD", "PNEBI_LUT"};

    AngularRadonSpectrum2d::AngularRadonSpectrum2d()
    : coeffs_(), isotropicKer_(), arsfOrder_(0),
    thetaToll_(M_PI / 180.0 * 0.5) { // pnebiLut_(), mode_(PNEBI_LUT), anisotropicStep_(720)
    }

    AngularRadonSpectrum2d::AngularRadonSpectrum2d(const std::vector<double>& coeffs)
    : coeffs_(coeffs), isotropicKer_(), arsfOrder_(0),
    thetaToll_(M_PI / 180.0 * 0.5) { // pnebiLut_(), mode_(PNEBI_LUT), anisotropicStep_(720)
    }

    AngularRadonSpectrum2d::~AngularRadonSpectrum2d() {
    }

    void AngularRadonSpectrum2d::setARSFOrder(int n) {
        arsfOrder_ = n;
        coeffs_.resize(2 * n + 2);
        //initLUT();
    }

    void AngularRadonSpectrum2d::setThetaToll(double thetaToll) {
        thetaToll_ = thetaToll;
    }

    void AngularRadonSpectrum2d::setComputeMode(ArsKernelIso2dComputeMode mode) {
        //mode_ = mode;
        isotropicKer_.setComputeMode(mode);
    }

    const std::string& AngularRadonSpectrum2d::getComputeModeName() const {
        //return MODE_NAME[mode_];
        return isotropicKer_.getComputeModeName();
    }

    const std::vector<double>& AngularRadonSpectrum2d::coefficients() const {
        return coeffs_;
    }

    double AngularRadonSpectrum2d::normCorr() const {
        double ret = 0.0;
        for (int i = 0; i < coeffs_.size(); ++i) {
            ret += coeffs_[i] * coeffs_[i];
        }
        return ret;
    }

    void AngularRadonSpectrum2d::setCoefficients(const std::vector<double>& coeffs) {
        coeffs_ = coeffs;
    }

    void AngularRadonSpectrum2d::setCoefficients(double* coeffs, size_t coeffsSz) {
        coeffs_.clear();
        for (size_t i = 0; i < coeffsSz; i++) {
            double c = coeffs[i];
            //            std::cout << "coeff[" << i << "] " << coeffs[i] << std::endl;
            coeffs_.push_back(c);
        }
    }

    void AngularRadonSpectrum2d::insertIsotropicGaussians(const VecVec2d& means, double sigma) {
        size_t kernelNum = means.size();
        double w = 1.0 / (kernelNum * kernelNum);
        //std::cout << "kernelNum " << kernelNum << ", mode_ " << mode_ << " " << MODE_NAME[mode_] << std::endl;

        //        if (pnebiLut_.getOrderMax() < arsfOrder_) {
        //            std::cerr << __FILE__ << "," << __LINE__ << ": LUT not initialized to right order. Initialized now." << std::endl;
        //            pnebiLut_.init(arsfOrder_, 0.0001);
        //        }

        if (coeffs_.size() != 2 * arsfOrder_ + 2) {
            coeffs_.resize(2 * arsfOrder_ + 2);
        }

        std::fill(coeffs_.begin(), coeffs_.end(), 0.0);

        for (int i = 0; i < kernelNum; ++i) {
            for (int j = i + 1; j < kernelNum; ++j) {
                cuars::ScopedTimer timer("ArsKernelIsotropic2d::computeFourier()");
                isotropicKer_.init(means[i], means[j], sigma);
                isotropicKer_.updateFourier(arsfOrder_, coeffs_, w);

            }
        }
    }

    void AngularRadonSpectrum2d::insertAnisotropicGaussians(const VecVec2d& means, const VecMat2d& covars, const std::vector<double>& weights) {
        ArsKernelAnisotropic2d nik;
        std::vector<double> coeffsPartial(arsfOrder_);
        int kernelNum = means.size();
        double wij;

        if (kernelNum != covars.size()) {
            std::cerr << __FILE__ << "," << __LINE__ << ": inconsistent vector sizes: found " << means.size()
                    << " mean values and " << covars.size() << " covariance matrices" << std::endl;
            return;
        }

        nik.setFourierOrder(arsfOrder_);
        //ARS_ASSERT(coeffs_.size() == 2 * arsfOrder_ && coeffsPartial.size() == 2 * arsfOrder_);
        coeffs_.resize(2 * arsfOrder_);
        coeffsPartial.resize(2 * arsfOrder_);


        std::fill(coeffs_.begin(), coeffs_.end(), 0.0);
        for (int i = 0; i < kernelNum; ++i) {
            for (int j = i + 1; j < kernelNum; ++j) {
                cuars::ScopedTimer timer("AnisotropicKernel::computeFourier()");
                nik.init(means[i], covars[i], means[j], covars[j]);
                nik.computeFourier(coeffsPartial);

                wij = weights[i] * weights[j];
                for (int f = 0; f < coeffs_.size(); ++f) {
                    coeffs_[f] += wij * coeffsPartial[f];
                }
            }
        }
    }

    void AngularRadonSpectrum2d::initLUT(double precision = 0.001) {
        //pnebiLut_.init(arsfOrder_, precision);
        isotropicKer_.initPnebiLut(arsfOrder_, precision);
    }

    double AngularRadonSpectrum2d::eval(double theta) const {
        return evaluateFourier(coeffs_, 2.0 * theta);
    }

    double AngularRadonSpectrum2d::findMax() const {
        double arsfMax, thetaMax;
        findGlobalMaxBBFourier(coeffs_, 0, M_PI, thetaToll_, 10.0, thetaMax, arsfMax);
        return arsfMax;
    }

    double AngularRadonSpectrum2d::findMax(double& thetaMax) const {
        double arsfMax;
        findGlobalMaxBBFourier(coeffs_, 0, M_PI, thetaToll_, 10.0, thetaMax, arsfMax);
        return arsfMax;
    }

    double AngularRadonSpectrum2d::findMax(double thetaLow, double thetaUpp) const {
        double arsfMax, thetaMax;
        findGlobalMaxBBFourier(coeffs_, thetaLow, thetaUpp, thetaToll_, 10.0, thetaMax, arsfMax);
        return arsfMax;
    }

    double AngularRadonSpectrum2d::findMax(double& thetaMax, double thetaLow, double thetaUpp) const {
        double arsfMax;
        findGlobalMaxBBFourier(coeffs_, thetaLow, thetaUpp, thetaToll_, 10.0, thetaMax, arsfMax);
        return arsfMax;
    }


} // end of namespace

