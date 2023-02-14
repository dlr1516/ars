#include <iostream>

#include <ars/BBTranslation.h>

int main (int argc, char **argv) {

    ars::VectorVector3 ptsA, ptsB;
    int dim = 3;
    Eigen::Affine2d transf;


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
    // ars::createRotationMatrix(transfTrue, 0, 0, 30 * M_PI / 180.0);
    transfTrue.translation() = Eigen::Vector3d(-1,-1,-1);
    std::cout << "Applying transformation transfTrue:" << std::endl
              << transfTrue.matrix() << std::endl;
    for (int i = 0; i < 6; ++i) //6 = number of points
    {
        Eigen::Vector3d ptsAiTransf = ptsA[i] + transfTrue.translation().matrix();
        ptsB.push_back(Eigen::Vector3d(ptsAiTransf(0), ptsAiTransf(1), ptsAiTransf(2)));
    }

    return 0;
}