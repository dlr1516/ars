#ifndef ARS_PROCRUSTES_UMEYAMA_H_
#define ARS_PROCRUSTES_UMEYAMA_H_

#include <vector>
#include <algorithm>

#include <ars/definitions.h>
#include <ars/utils.h>

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>

#define TWO_ 2
#define THREE_ 3

template <typename T>
using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename Scalar, int rank, typename sizeType>
auto Tensor_to_Matrix(const Eigen::Tensor<Scalar, rank> &tensor, const sizeType rows, const sizeType cols)
{
    return Eigen::Map<const MatrixType<Scalar>>(tensor.data(), rows, cols);
}

template <typename Scalar, typename... Dims>
auto Matrix_to_Tensor(const MatrixType<Scalar> &matrix, Dims... dims)
{
    constexpr int rank = sizeof...(Dims);
    return Eigen::TensorMap<Eigen::Tensor<const Scalar, rank>>(matrix.data(), {dims...});
}

inline void procrustesUmeyama3d(Eigen::Affine3d &transfOut, const ars::VectorVector3 &cloudA, const ars::VectorVector3 &cloudB)
{
    transfOut = Eigen::Affine3d::Identity();

    ARS_ASSERT(cloudA.size() == cloudB.size());

    // Input clouds/matrices are supposed to have size m x n
    int m = THREE_;
    int n = std::min<int>(cloudA.size(), cloudB.size()); // TODO: fix when size(cloudA)!=size(cloudB)

    // // cloudA->points.erase(cloudA->points.end() - 1);
    // // cloudB->points.erase(cloudB->points.end() - 1);
    // Eigen::MatrixXd clAMat = cloudA->getMatrixXdMap();
    // Eigen::MatrixXd clBMat = cloudB->getMatrixXdMap();
    // clAMat = clAMat.block<3,3>(0,0);
    // clBMat = clBMat.block<3,3>(0,0);
    Eigen::MatrixXd clAMat(m, n);
    Eigen::MatrixXd clBMat(m, n);

    for (int i = 0; i < n; ++i)
    {
        Eigen::Vector3d ptA = cloudA[i];
        Eigen::Vector3d ptB = cloudB[i];
        clAMat.col(i) = ptA;
        clBMat.col(i) = ptB;
    }

    std::cout << "m " << m << " n " << n << std::endl;
    std::cout << "Mat a rows " << clAMat.rows() << " cols " << clAMat.cols() << std::endl;
    std::cout << "clAMat" << std::endl
              << clAMat << std::endl;
    std::cout << "Mat b rows " << clBMat.rows() << " cols " << clBMat.cols() << std::endl;
    std::cout << "clBMat" << std::endl
              << clBMat << std::endl;

    Eigen::MatrixXd svd1InputMat = clAMat * clBMat.transpose();

    Eigen::JacobiSVD<Eigen::MatrixXd> svd1(svd1InputMat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd u1 = svd1.matrixU();
    Eigen::MatrixXd v1 = svd1.matrixV();

    Eigen::Matrix3d s_min = Eigen::Matrix3d::Identity();
    if (svd1InputMat.determinant() < 0)
        s_min(2, 2) = -1;

    // Note: for floating values matrices, you might get more accurate results with Eigen::ColPivHouseholderQR< MatrixType >
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> luDecomp(svd1InputMat); // needed to compute rank
    auto rank_abt = luDecomp.rank();
    if ((int)rank_abt < m - 1)
    {
        std::cerr << "Error: rank(a*b') < dim - 1 -> returning eye transf" << std::endl;
        // return;
    }

    Eigen::Matrix3d R1 = Eigen::Matrix3d::Identity();
    if (rank_abt == m - 1)
    {
        Eigen::Matrix3d s_argmin = Eigen::Matrix3d::Identity();
        if (u1.determinant() * v1.determinant() == -1)
            s_argmin(2, 2) = -1;
        R1 = u1 * s_argmin * v1.transpose();
    }
    else if (rank_abt == m)
        R1 = u1 * s_min * v1.transpose();

    // TODO: add computation of minimum value of mean squared error

    Eigen::Vector3d a_centroid = clAMat.rowwise().mean(); // mu_x
    Eigen::Vector3d b_centroid = clBMat.rowwise().mean(); // mu_y

    Eigen::MatrixXd a_centroid_repmat(a_centroid.rowwise().replicate(n));
    Eigen::MatrixXd diff_a = clAMat - a_centroid_repmat;
    Eigen::MatrixXd a_sqnorm = diff_a.colwise().squaredNorm();
    double sigma_a = a_sqnorm.sum() / n; // sigma_x

    Eigen::MatrixXd b_centroid_repmat(b_centroid.rowwise().replicate(n));
    Eigen::MatrixXd diff_b = clBMat - b_centroid_repmat;
    Eigen::MatrixXd b_sqnorm = diff_b.colwise().squaredNorm();
    double sigma_b = b_sqnorm.sum() / n; // sigma_y

    Eigen::Tensor<double, THREE_> sigma_xy_complete(m, m, n); // Sigma_xy
    sigma_xy_complete.setZero();
    for (int i = 0; i < n; ++i)
    {
        // MATLAB equivalent:
        // sigma_xy_complete.block<dim,dim,n>(0,0,i) = diff_b.col(i) * diff_a.col(ii).transpose();

        Eigen::Vector3d dAi = diff_a.col(i);
        Eigen::Vector3d dBi = diff_b.col(i);
        Eigen::MatrixXd tmp = dBi * dAi.transpose(); // this is actually a Matrix2d, but conversion function takes MatrixXd as input

        Eigen::Tensor<double, TWO_> tmpTens = Matrix_to_Tensor(tmp, THREE_, THREE_);
        sigma_xy_complete.chip(i, 2) = tmpTens;

        // std::cout << "i " << i << std::endl;
        // std::cout << "dAi " << std::endl
        //           << dAi << std::endl;
        // std::cout << "dBi " << std::endl
        //           << dBi << std::endl;
        // std::cout << "tmp " << std::endl
        //           << tmp << std::endl;
        // std::cout << "sigma_xy_complete.chip(i, 2) " << std::endl
        //           << sigma_xy_complete.chip(i, 2) << std::endl;
        // std::cout << "sigma_xy_complete " << std::endl
        //           << sigma_xy_complete << std::endl;
        // std::cout << "tmpTens " << std::endl
        //           << tmpTens << std::endl;
        // std::cout << std::endl;
    }
    // std::cout << "sigma_xy_complete " << std::endl
    //           << sigma_xy_complete << std::endl;
    Eigen::Tensor<double, 2> sigma_xy_tensor(THREE_, THREE_);
    sigma_xy_tensor.setZero(); //!!
    // built-in version
    // std::cout << std::endl << "sigma_xy_tensor" << sigma_xy_tensor << std::endl;
    // std::array<int, 3> two_dims{{3,3}};
    // sigma_xy_tensor = sigma_xy_complete.sum(two_dims);
    // hand-made version
    for (int i = 0; i < n; ++i)
    {
        sigma_xy_tensor += sigma_xy_complete.chip(i, 2); // this still needs normalization!
    }
    // std::cout << "sigma_xy_tensor " << std::endl
    //           << sigma_xy_tensor << std::endl; // this still needs normalization!
    Eigen::MatrixXd sigma_xy = Tensor_to_Matrix(sigma_xy_tensor, m, m) / n;
    // std::cout << "sigma_xy " << std::endl
    //   << sigma_xy << std::endl;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd_sigma_xy(sigma_xy, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::FullPivLU<Eigen::Matrix3d> sigma_xy_decomp(sigma_xy); // needed to compute rank
    auto rank_sigma_xy = sigma_xy_decomp.rank();
    if ((int)rank_sigma_xy < m - 1)
    {
        std::cerr << "Error: rank(sigma_xy) < dim - 1 -> returning eye() transf" << std::endl;
        return;
    }

    Eigen::MatrixXd u_sigma_xy = svd_sigma_xy.matrixU();
    // Eigen::DiagonalMatrix<double, TWO_>(svd_sigma_xy.singularValues())
    Eigen::MatrixXd d_sigma_xy = svd_sigma_xy.singularValues().asDiagonal(); // called S in Eigen documentation
    Eigen::MatrixXd v_sigma_xy = svd_sigma_xy.matrixV();

    Eigen::Matrix3d s_sigma_xy = Eigen::Matrix3d::Identity();
    if (sigma_xy.determinant() < 0)
        s_sigma_xy(1, 1) = -1;
    if ((int)rank_sigma_xy == m - 1)
        if (u_sigma_xy.determinant() * v_sigma_xy.determinant() == -1)
            s_sigma_xy(1, 1) = -1;

    Eigen::Matrix3d R2 = u_sigma_xy * s_sigma_xy * v_sigma_xy.transpose();
    double c = (d_sigma_xy * s_sigma_xy).trace() / sigma_a;
    Eigen::Vector3d transl = b_centroid - c * R2 * a_centroid;

    transfOut.linear() = R2;
    transfOut.translation() = transl;
    transfOut.makeAffine();
}

inline void procrustesUmeyama2d(Eigen::Affine2d &transfOut, const ars::VectorVector2 &cloudA, const ars::VectorVector2 &cloudB)
{
    transfOut = Eigen::Affine2d::Identity();

    ARS_ASSERT(cloudA.size() == cloudB.size());

    // Input clouds/matrices are supposed to have size m x n
    int m = TWO_;
    int n = std::min<int>(cloudA.size(), cloudB.size()); // TODO: fix when size(cloudA)!=size(cloudB)

    Eigen::MatrixXd clAMat(m, n);
    Eigen::MatrixXd clBMat(m, n);

    for (int i = 0; i < n; ++i)
    {
        Eigen::Vector2d ptA = cloudA[i];
        Eigen::Vector2d ptB = cloudB[i];
        clAMat.col(i) = ptA;
        clBMat.col(i) = ptB;
    }

    // std::cout << "m " << m << " n " << n << std::endl;
    // std::cout << "Mat a rows " << clAMat.rows() << " cols " << clAMat.cols() << std::endl;
    // std::cout << "clAMat" << std::endl
    //           << clAMat << std::endl;
    // std::cout << "Mat b rows " << clBMat.rows() << " cols " << clBMat.cols() << std::endl;
    // std::cout << "clBMat" << std::endl
    //           << clBMat << std::endl;

    Eigen::MatrixXd svd1InputMat = clAMat * clBMat.transpose();

    Eigen::JacobiSVD<Eigen::MatrixXd> svd1(svd1InputMat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd u1 = svd1.matrixU();
    Eigen::MatrixXd v1 = svd1.matrixV();

    Eigen::Matrix2d s_min = Eigen::Matrix2d::Identity();
    if (svd1InputMat.determinant() < 0)
        s_min(1, 1) = -1;

    // Note: for floating values matrices, you might get more accurate results with Eigen::ColPivHouseholderQR< MatrixType >
    Eigen::FullPivLU<Eigen::MatrixXd> luDecomp(svd1InputMat); // needed to compute rank
    auto rank_abt = luDecomp.rank();
    if ((int)rank_abt < m - 1)
    {
        std::cerr << "Error: rank(a*b') < dim - 1 -> returning eye transf" << std::endl;
        return;
    }

    Eigen::Matrix2d R1 = Eigen::Matrix2d::Identity();
    if (rank_abt == m - 1)
    {
        Eigen::Matrix2d s_argmin = Eigen::Matrix2d::Identity();
        if (u1.determinant() * v1.determinant() == -1)
            s_argmin(1, 1) = -1;
        R1 = u1 * s_argmin * v1.transpose();
    }
    else if (rank_abt == m)
        R1 = u1 * s_min * v1.transpose();

    // TODO: add computation of minimum value of mean squared error

    // std::cout << "m " << m << " n " << n << std::endl;
    // std::cout << "Mat a rows " << clAMat.rows() << " cols " << clAMat.cols() << std::endl;

    Eigen::Vector2d a_centroid = clAMat.rowwise().mean(); // mu_x
    Eigen::Vector2d b_centroid = clBMat.rowwise().mean(); // mu_y
    // std::cout << "a_centroid " << std::endl
    //           << a_centroid << std::endl;

    Eigen::MatrixXd a_centroid_repmat(a_centroid.rowwise().replicate(n));
    // std::cout << "a_centroid_repmat rows " << a_centroid_repmat.rows() << " cols " << a_centroid_repmat.cols() << std::endl;
    // std::cout << std::endl << "a_centroid_repmat " << std::endl << a_centroid_repmat << std::endl;
    Eigen::MatrixXd diff_a = clAMat - a_centroid_repmat;
    // std::cout << std::endl << "diff_a " << std::endl << diff_a << std::endl;
    Eigen::MatrixXd a_sqnorm = diff_a.colwise().squaredNorm();
    // std::cout << std::endl << "a_sqnorm " << std::endl << a_sqnorm << std::endl;
    double sigma_a = a_sqnorm.sum() / n; // sigma_x
    // std::cout << "sigma_a " << sigma_a << std::endl;

    Eigen::MatrixXd b_centroid_repmat(b_centroid.rowwise().replicate(n));
    // std::cout << "b_centroid_repmat rows " << b_centroid_repmat.rows() << " cols " << b_centroid_repmat.cols() << std::endl;
    // std::cout << std::endl
    //           << "b_centroid_repmat " << std::endl
    //           << b_centroid_repmat << std::endl;
    Eigen::MatrixXd diff_b = clBMat - b_centroid_repmat;
    // std::cout << std::endl
    //           << "diff_b " << std::endl
    //           << diff_b << std::endl;
    Eigen::MatrixXd b_sqnorm = diff_b.colwise().squaredNorm();
    // std::cout << std::endl
    //           << "b_sqnorm " << std::endl
    //           << b_sqnorm << std::endl;
    double sigma_b = b_sqnorm.sum() / n; // sigma_y
    // std::cout << "sigma_b " << sigma_b << std::endl;

    Eigen::Tensor<double, THREE_> sigma_xy_complete(m, m, n); // Sigma_xy
    sigma_xy_complete.setZero();
    for (int i = 0; i < n; ++i)
    {
        // MATLAB equivalent:
        // sigma_xy_complete.block<dim,dim,n>(0,0,i) = diff_b.col(i) * diff_a.col(ii).transpose();

        Eigen::Vector2d dAi = diff_a.col(i);
        Eigen::Vector2d dBi = diff_b.col(i);
        Eigen::MatrixXd tmp = dBi * dAi.transpose(); // this is actually a Matrix2d, but conversion function takes MatrixXd as input

        Eigen::Tensor<double, TWO_> tmpTens = Matrix_to_Tensor(tmp, TWO_, TWO_);
        sigma_xy_complete.chip(i, 2) = tmpTens;

        // std::cout << "i " << i << std::endl;
        // std::cout << "dAi " << std::endl
        //           << dAi << std::endl;
        // std::cout << "dBi " << std::endl
        //           << dBi << std::endl;
        // std::cout << "tmp " << std::endl
        //           << tmp << std::endl;
        // std::cout << "sigma_xy_complete.chip(i, 2) " << std::endl
        //           << sigma_xy_complete.chip(i, 2) << std::endl;
        // std::cout << "sigma_xy_complete " << std::endl
        //           << sigma_xy_complete << std::endl;
        // std::cout << "tmpTens " << std::endl
        //           << tmpTens << std::endl;
        // std::cout << std::endl;
    }
    // std::cout << "sigma_xy_complete " << std::endl
    //           << sigma_xy_complete << std::endl;
    Eigen::Tensor<double, 2> sigma_xy_tensor(TWO_, TWO_);
    sigma_xy_tensor.setZero(); //!!
    // built-in version
    // std::cout << std::endl << "sigma_xy_tensor" << sigma_xy_tensor << std::endl;
    // std::array<int, 3> two_dims{{2,2}};
    // sigma_xy_tensor = sigma_xy_complete.sum(two_dims);
    // hand-made version
    for (int i = 0; i < n; ++i)
    {
        sigma_xy_tensor += sigma_xy_complete.chip(i, 2); // this still needs normalization!
    }
    // std::cout << "sigma_xy_tensor " << std::endl
    //   << sigma_xy_tensor << std::endl; // this still needs normalization!
    Eigen::MatrixXd sigma_xy = Tensor_to_Matrix(sigma_xy_tensor, m, m) / n;
    // std::cout << "sigma_xy " << std::endl
    //   << sigma_xy << std::endl;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd_sigma_xy(sigma_xy, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::FullPivLU<Eigen::Matrix2d> sigma_xy_decomp(sigma_xy); // needed to compute rank
    auto rank_sigma_xy = sigma_xy_decomp.rank();
    if ((int)rank_sigma_xy < m - 1)
    {
        std::cerr << "Error: rank(sigma_xy) < dim - 1 -> returning eye() transf" << std::endl;
        return;
    }

    Eigen::MatrixXd u_sigma_xy = svd_sigma_xy.matrixU();
    // Eigen::DiagonalMatrix<double, TWO_>(svd_sigma_xy.singularValues())
    Eigen::MatrixXd d_sigma_xy = svd_sigma_xy.singularValues().asDiagonal(); // called S in Eigen documentation
    Eigen::MatrixXd v_sigma_xy = svd_sigma_xy.matrixV();

    Eigen::Matrix2d s_sigma_xy = Eigen::Matrix2d::Identity();
    if (sigma_xy.determinant() < 0)
        s_sigma_xy(1, 1) = -1;
    if ((int)rank_sigma_xy == m - 1)
        if (u_sigma_xy.determinant() * v_sigma_xy.determinant() == -1)
            s_sigma_xy(1, 1) = -1;

    Eigen::Matrix2d R2 = u_sigma_xy * s_sigma_xy * v_sigma_xy.transpose();
    double c = (d_sigma_xy * s_sigma_xy).trace() / sigma_a;
    Eigen::Vector2d transl = b_centroid - c * R2 * a_centroid;

    transfOut.linear() = R2;
    transfOut.translation() = transl;
    transfOut.makeAffine();
}

#endif /*ARS_PROCRUSTES_UMEYAMA_H_*/