#include "utils.h"
#include "definitions.h"

#define INV_SQRT_2_PI 0.3989422804 // 1 / sqrt(2.0 * M_PI)

struct ArsKernelAnisotropic2d_simpl {
    int nFourier_;
    double muMod_;
    double muAng_;
    double sigmaMod_;
    double sigmaAng_;
    double sigmaDif_;

    std::vector<double> kernelVal_;
    std::vector<std::complex<double> > freqvec_;

    struct DataLut {
        double varCos;
        double varSin;
        double meanConst;
        double meanCos;
        double meanSin;
        std::vector<double> cosTh;
        std::vector<double> sinTh;
    };

    DataLut lut_;

    void init(const cuars::Vec2d &mean1, const cuars::Mat2d &covar1, const cuars::Vec2d &mean2, const cuars::Mat2d &covar2) {
        cuars::Vec2d mu12;
        cuars::Mat2d sigma12;
        double a, b, lmax, lmin, c, s;

        cuars::vec2diff(mu12, mean2, mean1);
        //        muMod_ = mu12.norm();
        muMod_ = cuars::vec2norm(mu12);
        //        muAng_ = atan2(mu12.data_[1], mu12.data_[0]);
        muAng_ = atan2(mu12.y, mu12.x);

        cuars::mat2dSum(sigma12, covar1, covar2);

        // Diagonalizes sigma12
        cuars::diagonalize(sigma12, lmin, lmax, sigmaAng_);

        //        a = 0.5 * (sigma12(1, 1) - sigma12(0, 0));
        //        b = 0.5 * (sigma12(0, 1) + sigma12(1, 0));
        //        ARS_VARIABLE2(a, b);
        //
        //        sigmaAng_ = 0.5 * atan2(-b, a);
        //
        //        c = cos(sigmaAng_);
        //        s = sin(sigmaAng_);
        //        lmax = sigma12(0, 0) * c * c + sigma12(1, 1) * s * s + (sigma12(0, 1) + sigma12(1, 0)) * c * s;
        //        lmin = sigma12(0, 0) * s * s + sigma12(1, 1) * c * c - (sigma12(0, 1) + sigma12(1, 0)) * c * s;
        //        ARS_VARIABLE3(sigmaAng_, lmax, lmin);
        //
        //        if (lmax < lmin) {
        //            sigmaAng_ += 0.5 * M_PI;
        //            std::swap(lmax, lmin);
        //            ARS_PRINT("lmin " << lmin << " < lmax " << lmax << ": swap, sigmaAng_ + PI/2: " << sigmaAng_);
        //        }

        sigmaMod_ = 0.5 * (lmax + lmin);
        sigmaDif_ = (lmax - lmin) / (lmax + lmin);

        //        ARS_PRINT("muMod_ " << muMod_ << ", muAng_[rad] " << muAng_ << " [deg] " << (180.0 / M_PI * muAng_) << "\n"
        //                << "sigmaMod_ " << sigmaMod_ << ", sigmaAng_[rad] " << sigmaAng_ << " [deg] " << (180.0 / M_PI * sigmaAng_)
        //                << ", sigmaDif_ " << sigmaDif_ << "\n");

        // Initizalizes LUT
        lut_.varCos = sigmaMod_ * sigmaDif_ * cos(2.0 * sigmaAng_);
        lut_.varSin = sigmaMod_ * sigmaDif_ * sin(2.0 * sigmaAng_);
        lut_.meanConst = 0.5 * muMod_ * muMod_;
        lut_.meanCos = lut_.meanConst * cos(2.0 * muAng_);
        lut_.meanSin = lut_.meanConst * sin(2.0 * muAng_);
    }

    void initCosSinLut() {
        double dt = M_PI / nFourier_;
        int n2 = nFourier_ << 1;
        lut_.cosTh.resize(n2);
        lut_.sinTh.resize(n2);
        for (int i = 0; i < n2; ++i) {
            lut_.cosTh[i] = cos(dt * i); // cosTh[i] = cos(2 * PI * i / (2 * nFourier)) = cos(PI / nFourier * i)
            lut_.sinTh[i] = sin(dt * i); // sinTh[i] = sin(2 * PI * i / (2 * nFourier)) = sin(PI / nFourier * i)
        }
        // kernelVal_ stores samples of the ARS kernel sampled at frequency double than DFT/FFT order
        kernelVal_.resize(n2);
        freqvec_.resize(nFourier_ + 1);
    }

    void setFourierOrder(int nFourier) {
        if (nFourier_ != nFourier) {
            nFourier_ = nFourier;
            initCosSinLut();
        }
    }

    Eigen::FFT<double> fft_;

    void computeFourier(std::vector<double> &coeffs) {
        //double sumCos, sumSin, cosCurr, cosNext, cosIncr, sinCurr, sinNext, sinIncr;
        double dt, factor, varInv;
        int n2 = nFourier_ << 1;

        ARS_ASSERT(nFourier_ > 0);

        // Evaluates each of the Fourier coefficients
        if (coeffs.size() != n2) {
            coeffs.resize(n2);
        }

        // Evaluates the kernel function at given intervals
        //        dt = M_PI / nFourier_;
        //        for (int i = 0; i < nFourier; ++i) {
        //            kernelVal[i] = value(dt * i);
        //        }
        // Alternative and faster evaluation using LUT
        //  var = sigmaMod_ * (1.0 + sigmaDif_ * cos(2.0 * t - 2.0 * sigmaAng_))
        //      = sigmaMod_ + sigmaMod_ * sigmaDif_ * (cos(2.0 * t) * cos(2.0 * sigmaAng_) + sin(2.0 * t) * sin(2.0 * sigmaAng_))
        //      = sigmaMod_ + (sigmaMod_ * sigmaDif_ * cos(2.0 * sigmaAng_)) * cos(2.0 * t) + (sigmaMod_ * sigmaDif_ * sin(2.0 * sigmaAng_)) * sin(2.0 * t)
        //      = sigmaMod_ + lut_.varCos * cosTh[i] + lut_.varSin * lut_.sinTh[i]
        for (int i = 0; i < n2; ++i) {
            varInv = 1.0 / (sigmaMod_ + lut_.varCos * lut_.cosTh[i] + lut_.varSin * lut_.sinTh[i]);
            kernelVal_[i] = INV_SQRT_2_PI * sqrt(varInv) * exp(-0.5 * (lut_.meanConst + lut_.meanCos * lut_.cosTh[i] + lut_.meanSin * lut_.sinTh[i]) * varInv);
        }

        //fft(kernelVal, coeffs, nFourier);

        fft_.fwd(freqvec_, kernelVal_);

        coeffs[0] = 0.5 * freqvec_[0].real() / nFourier_;
        coeffs[1] = 0.0;
        factor = 1.0 / nFourier_;
        for (int i = 1; i < nFourier_; ++i) {
            coeffs[2 * i] = factor * freqvec_[i].real();
            coeffs[2 * i + 1] = -factor * freqvec_[i].imag();
        }
    }
};

__global__
void ianigK() {
    printf("Hello from InsertAnisotropicGaussians kernel\n");
}

__host__
void insertAnisotropicGaussians_Clang(int arsfOrder_, const cuars::VecVec2d& means, const cuars::VecMat2d& covars, const std::vector<double>& weights, std::vector<double>& coeffs_) {
    ArsKernelAnisotropic2d_simpl nik;
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
            //            cuars::ScopedTimer timer("AnisotropicKernel::computeFourier()");
            nik.init(means[i], covars[i], means[j], covars[j]);
            nik.computeFourier(coeffsPartial);

            wij = weights[i] * weights[j];
            for (int f = 0; f < coeffs_.size(); ++f) {
                coeffs_[f] += wij * coeffsPartial[f];
            }
        }
    }

}
