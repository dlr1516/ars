#include <iostream>

// #include ""

#include "ars/ProcrustesUmeyama.h"

// % 3D setup
// a = [1 1 1; 3 4 5; 10 8 6; 2 4 7; 7 1 2; -14 -5 16]';
// b = zeros(size(a));
// rotz_true = 30; % deg
// transl_true = [-1; -1; -1];
// for ii = 1:size(a, 2)
//     b(:, ii) = rotz(rotz_true) * a(:, ii) + transl_true;
// end
// disp(b);
// d = 3;

// %umeyama paper setup
// % a = [0 0; 1 0 ; 0 2]';
// % b = [0 0; -1 0; 0 2]';
// % d = 2;

// rigid_out = procrustes_umeyama(a,b,d);
// disp("transf_out");
// disp(rigid_out.A)

void create_rotation_matrix(Eigen::Affine3d &rotM, double ax, double ay, double az)
{
    Eigen::Affine3d rx =
        Eigen::Affine3d(Eigen::AngleAxisd(ax, Eigen::Vector3d(1, 0, 0)));
    Eigen::Affine3d ry =
        Eigen::Affine3d(Eigen::AngleAxisd(ay, Eigen::Vector3d(0, 1, 0)));
    Eigen::Affine3d rz =
        Eigen::Affine3d(Eigen::AngleAxisd(az, Eigen::Vector3d(0, 0, 1)));
    rotM = rz * ry * rx;
}

int main(int argc, char **argv)
{

    // START OF 3D TEST
    ars::VectorVector3 ptsA, ptsB;
    int dim = 3;
    Eigen::Affine3d transf;

    // for Umeyama paper setup:
    // Eigen::Vector3d aA(Eigen::Vector3d(0,0,0));
    // ptsA.push_back(aA);
    // Eigen::Vector3d bA(Eigen::Vector3d(1,0,0));
    // ptsA->push_back(bA);
    // Eigen::Vector3d cA(Eigen::Vector3d(0,2,0));
    // ptsA->push_back(cA);
    // Eigen::Vector3d aB(Eigen::Vector3d(0,0,0));
    // ptsB->push_back(aB);
    // Eigen::Vector3d bB(Eigen::Vector3d(-1,0,0));
    // ptsB->push_back(bB);
    // Eigen::Vector3d cB(Eigen::Vector3d(0,2,0));
    // ptsB->push_back(cB);

    // test with other points:
    Eigen::Vector3d aA(Eigen::Vector3d(1, 1, 1));
    Eigen::Vector3d bA(Eigen::Vector3d(3, 4, 5));
    Eigen::Vector3d cA(Eigen::Vector3d(10, 8, 6));
    Eigen::Vector3d dA(Eigen::Vector3d(2, 4, 7));
    Eigen::Vector3d eA(Eigen::Vector3d(7, 1, 2));
    Eigen::Vector3d fA(Eigen::Vector3d(-14, -5, 16));
    Eigen::Vector3d zeroB(Eigen::Vector3d(0, 0, 0));
    ptsA.push_back(aA);
    ptsA.push_back(bA);
    ptsA.push_back(cA);
    ptsA.push_back(dA);
    ptsA.push_back(eA);
    ptsA.push_back(fA);
    Eigen::Affine3d transfTrue;
    create_rotation_matrix(transfTrue, 0, 0, 30 * M_PI / 180.0);
    transfTrue.translation() = Eigen::Vector3d(-1,-1,-1);
    std::cout << "Applying transformation transfTrue:" << std::endl
              << transfTrue.matrix() << std::endl;
    for (int i = 0; i < 6; ++i) //6 = number of points
    {
        Eigen::Vector3d ptsAiTransf = transfTrue.linear().matrix() * ptsA[i] + transfTrue.translation().matrix();
        ptsB.push_back(Eigen::Vector3d(ptsAiTransf(0), ptsAiTransf(1), ptsAiTransf(2)));
    }
    procrustesUmeyama3d(transf, ptsA, ptsB);
    // END OF 3D TEST

    // START OF 2D TEST
    // pcl::PointCloud<pcl::PointXY>::Ptr ptsA(new pcl::PointCloud<pcl::PointXY>), ptsB(new pcl::PointCloud<pcl::PointXY>);
    // int dim = 2;
    // Eigen::Affine2d transf;
    // pcl::PointXY aA;
    // aA.x = 0.0; aA.y = 0.0;
    // ptsA->push_back(aA);
    // pcl::PointXY bA;
    // bA.x = 1.0; bA.y = 0.0;
    // ptsA->push_back(bA);
    // pcl::PointXY cA;
    // cA.x = 0.0; cA.y = 2.0;
    // ptsA->push_back(cA);
    // pcl::PointXY aB;
    // aB.x = 0.0; aB.y = 0.0;
    // ptsB->push_back(aB);
    // pcl::PointXY bB;
    // bB.x = -1.0; bB.y = 0.0;
    // ptsB->push_back(bB);
    // pcl::PointXY cB;
    // cB.x = 0.0; cB.y = 2.0;
    // ptsB->push_back(cB);
    // procrustesUmeyama2d(transf, ptsA, ptsB);
    // END OF 2D TEST

    std::cout << "transfOut" << std::endl
              << transf.matrix() << std::endl;

    return 0;
}