#ifndef ARS_PROCRUSTES_UMEYAMA_H_
#define ARS_PROCRUSTES_UMEYAMA_H_

#include <vector>
#include <algorithm>

#include <ars/definitions.h>

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>

#define TWO_ 2
#define THREE_ 3

template<typename T>
using  MatrixType = Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>;
   
template<typename Scalar,int rank, typename sizeType>
auto Tensor_to_Matrix(const Eigen::Tensor<Scalar,rank> &tensor,const sizeType rows,const sizeType cols)
{
    return Eigen::Map<const MatrixType<Scalar>> (tensor.data(), rows,cols);
}


template<typename Scalar, typename... Dims>
auto Matrix_to_Tensor(const MatrixType<Scalar> &matrix, Dims... dims)
{
    constexpr int rank = sizeof... (Dims);
    return Eigen::TensorMap<Eigen::Tensor<const Scalar, rank>>(matrix.data(), {dims...});
}


void procrustesUmeyama3d(Eigen::Affine3d& transfOut, const ars::VectorVector3& cloudA, const ars::VectorVector3& cloudB);

void procrustesUmeyama2d(Eigen::Affine2d& transfOut, const ars::VectorVector2& cloudA, const ars::VectorVector2& cloudB);


// function transf_out = procrustes_umeyama(a,b, dim)
//     % Least squares estimation of transformation parameters between two
//     % point patterns (clouds) a, b with dimensionality dim = 2 or 3
//     % following procedure described in T-PAMI article by Shinji Umeyama, 1991
    
//     m = dim; % = size(a, 1)
//     n = size(a,2); % = number of points

//     if dim == 2
//         transf_out = rigidtform2d(0, [0,0]);
//     elseif dim == 3
//         transf_out = rigidtform3d([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
//     else
//         disp("Error: dim has to be 2 or 3!");
//         return;
//     end
    

//     %LEMMA
//     %     [u,d,v] = svd(a * b');
//     %     disp ("u")
//     %     disp (u)
//     %     disp ("d") %should be granted that diagonal elements are in decreasing order
//     %     disp (d)
//     %     disp ("v")
//     %     disp (v) 
//     %     
//     %     s_min = eye(dim);
//     %     if det(a * b') < 0 
//     %         s_min(dim,dim) = -1;
//     %     end
//     %     
//     %     rank_abt = rank(a*b');
//     %     if rank_abt < m - 1
//     %         disp("Error: rank(a*b') < dim - 1 -> returning eye transf");
//     %         return;
//     %     end
//     % 
//     %     R1 = eye(dim);
//     % 
//     %     if rank_abt == m-1
//     %         s_argmin = eye(dim);
//     %         if det(u) * det(v) == -1 
//     %             s_argmin(dim,dim) = -1;
//     %         end
//     %         R1 = u * s_argmin * v';
//     %     else %if rank_abt == m
//     %         R1 = u * s_min * v';
//     %     end
//     % 
//     %     disp("R1");
//     %     disp(R1);

//     %TODO: add computation of minimum value of mean squared error

//     a_centroid = mean(a,2); %mu_x
//     b_centroid = mean(b,2); %mu_y

//     a_centroid_vec = repmat(a_centroid,1,n);
//     diff_a = a - a_centroid_vec;
//     a_sqnorm = (vecnorm(diff_a,2,1)).^2;    
//     sigma_a = sum(a_sqnorm) / n; %sigma_x

//     b_centroid_vec = repmat(b_centroid,1,n);
//     diff_b = b - b_centroid_vec;
//     b_sqnorm = (vecnorm(diff_b,2,1)).^2;    
//     sigma_b = sum(b_sqnorm) / n; %sigma_y

//     sigma_xy_complete = zeros(dim,dim,n); %Sigma_xy
//     for ii = 1:n
//         sigma_xy_complete(:,:,ii) = diff_b(:,ii) * diff_a(:,ii)';
//     end
//     sigma_xy = sum(sigma_xy_complete, 3) ./ n;

//     [u_sigma_xy,d_sigma_xy,v_sigma_xy] = svd(sigma_xy);

//     rank_sigma_xy = rank(sigma_xy)
//     if rank_sigma_xy < m - 1
//         disp("Error: rank(sigma_xy) < dim - 1 -> returning eye() transf");
//         return;
//     end

//     s_sigma_xy = eye(dim);
//     if det(sigma_xy) < 0
//         s_sigma_xy (dim,dim) = -1;
//     end    
//     if rank_sigma_xy == m-1
//         if det(u_sigma_xy) * det(v_sigma_xy) == -1 
//             s_sigma_xy(dim,dim) = -1;
//         end
//     end
    

//     R2 = u_sigma_xy * s_sigma_xy * v_sigma_xy'
//     c = trace(d_sigma_xy * s_sigma_xy) / sigma_a
//     transl = b_centroid - c * R2 * a_centroid
    
//     transf_out.R = R2;
//     transf_out.Translation = transl;

// end

#endif /*CUDARS_PROCRUSTES_UMEYAMA_H_*/