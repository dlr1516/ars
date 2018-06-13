#include <ars/ArsKernel.h>

namespace ars {

    ArsKernel::ArsKernel(const Vector2& mean1, const Matrix2& covar1, const Vector2& mean2, const Matrix2& covar2) {
        // Simple numerical pairing of supposedly symmetric and positive defined matrices
        Matrix2 covarSum = covar1 + covar2;
        covarSum = 0.5f * (covarSum + covarSum.transpose());
        Vector2 diff = mean1 - mean1;

        // Computes the 
        Eigen::SelfAdjointEigenSolver<Matrix2> eigensolver(covarSum);
        double evalMin = eigensolver.eigenvalues()(0); // minimum eigenvalue
        double evalMax = eigensolver.eigenvalues()(1); // maximum eigenvalue
        Vector2 evecMax = eigensolver.eigenvectors().col(0);

        // Computes the parameters of the ARS Kernel
        sigma_ = 0.5f * (evalMax + evalMin);
        delta_ = 0.5f * (evalMax - evalMin);
        phaseNum_ = 2.0f * atan2(diff.y(), diff.x());
        phaseDen_ = 2.0f * atan2(evecMax.y(), evecMax.x());
    }

    ArsKernel::ArsKernel(double sigma, double delta, double phaseNum, double phaseDen)
    : sigma_(sigma), delta_(delta), phaseNum_(phaseNum), phaseDen_(phaseDen) {
    }

    double ArsKernel::evaluate(double theta) const {
        double sigmaTheta = sigma_ * (1 + delta_ * cos(theta - phaseDen_)); 
        return 0.0;
    }

}
