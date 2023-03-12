#include <iostream>

#include <ars/ProcrustesUmeyama.h>

// % 3D setup - MATLAB CODE
// a = [1 1 1; 3 4 5; 10 8 6; 2 4 7; 7 1 2; -14 -5 16]';
// b = zeros(size(a));
// rotz_true = 30; % deg
// transl_true = [-1; -1; -1];
// for ii = 1:size(a, 2)
//     b(:, ii) = rotz(rotz_true) * a(:, ii) + transl_true;
// end
// disp(b);
// d = 3;
// rigid_out = procrustes_umeyama(a,b,d);
// disp("transf_out");
// disp(rigid_out.A)


int main(int argc, char **argv)
{
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
    ars::createRotationMatrix(transfTrue, 0, 0, 30 * M_PI / 180.0);
    transfTrue.translation() = Eigen::Vector3d(-1,-1,-1);
    std::cout << "Applying transformation transfTrue:" << std::endl
              << transfTrue.matrix() << std::endl;
    for (int i = 0; i < 6; ++i) //6 = number of points
    {
        Eigen::Vector3d ptsAiTransf = transfTrue.linear().matrix() * ptsA[i] + transfTrue.translation().matrix();
        ptsB.push_back(Eigen::Vector3d(ptsAiTransf(0), ptsAiTransf(1), ptsAiTransf(2)));
    }
    procrustesUmeyama3d(transf, ptsA, ptsB);
        

    std::cout << "transfOut" << std::endl
              << transf.matrix() << std::endl;

    return 0;
}