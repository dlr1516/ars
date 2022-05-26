#include <ars/TranslationRefiner.h>

namespace ars {

    void TranslationRefiner::associate() {
        std::cout << "" << std::endl;
    }

    void TranslationRefiner::computeProcrustes(const ars::VectorVector2 &pointsSrc, const ars::VectorVector2 &pointsDst, Eigen::Affine2d &transf) {
        // Procustes
        ars::Vector2 meanSrc, meanDst;
        Eigen::MatrixXd S = Eigen::MatrixXd::Zero(2, 2);
        int n = std::min(pointsSrc.size(), pointsDst.size()); // works for now, counting on dummy examples not to have fake correspondences
        int numAssociations = 0;
        int sizeDst = pointsDst.size();
        meanSrc = ars::Vector2::Zero(2, 1);
        meanDst = ars::Vector2::Zero(2, 1);
        for (int i = 0; i < n; ++i) {
            // if (pointsDst[i].associated) {
            meanSrc += pointsSrc[i];
            meanDst += pointsDst[i];
            numAssociations++;
            // }
        }
        std::cout << "numAssociations " << numAssociations << std::endl;
        if (numAssociations != 0) {
            meanSrc = meanSrc * 1.0f / numAssociations;
            meanDst = meanDst * 1.0f / numAssociations;
        } else {
            return;
        }

        for (int i = 0; i < n && i < sizeDst; ++i) {
            // if (!pointsDst[i].associated)
            // continue;
            S += (pointsSrc[i] - meanSrc) * (pointsDst[i] - meanDst).transpose();
        }
        std::cout << "S before SVD" << std::endl
                << S << std::endl;

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(S, Eigen::ComputeThinU | Eigen::ComputeThinV);
        float d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
        if (d > 0)
            d = 1.0;
        else
            d = -1.0;
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(2, 2);
        I(1, 1) = d;
        Eigen::MatrixXd R = svd.matrixV() * I * svd.matrixU().transpose();

        std::cout << "S after SVD" << std::endl
                << S << std::endl;
        std::cout << "U after SVD" << std::endl
                << S << std::endl;
        std::cout << "V after SVD" << std::endl
                << S << std::endl;

        ars::Vector2 transl = ars::Vector2::Zero(2, 1);
        transl += meanDst - R * meanSrc;

        transf.linear() = R;
        transf.translation() = (meanDst - R * meanSrc);
        transf.makeAffine();
    }

    void TranslationRefiner::icp() {
        std::cout << "" << std::endl;
    }

    void TranslationRefiner::setMaxIterations(int numMaxIter) {
        maxIterations_ = numMaxIter;
    }

    void TranslationRefiner::setStopThresh(double th) {
        stopThresh_ = th;
    }



}