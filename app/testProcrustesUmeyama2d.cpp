#include <iostream>

#include <ars/ProcrustesUmeyama.h>

// %umeyama paper setup
// % a = [0 0; 1 0 ; 0 2]';
// % b = [0 0; -1 0; 0 2]';
// % d = 2;

int main(int argc, char **argv)
{
    ars::VectorVector2 ptsA, ptsB;
    int dim = 2;
    Eigen::Affine2d transf;

    // insert points for Umeyama paper setup
    Eigen::Vector2d aA(Eigen::Vector2d(0.0, 0.0));
    ptsA.push_back(aA);
    Eigen::Vector2d bA(Eigen::Vector2d(1.0, 0.0));
    ptsA.push_back(bA);
    Eigen::Vector2d cA(Eigen::Vector2d(0.0, 2.0));
    ptsA.push_back(cA);
    Eigen::Vector2d aB(Eigen::Vector2d(0.0, 0.0));
    ptsB.push_back(aB);
    Eigen::Vector2d bB(Eigen::Vector2d(-1.0,0.0));
    ptsB.push_back(bB);
    Eigen::Vector2d cB(Eigen::Vector2d(0.0,2.0));
    ptsB.push_back(cB);
    procrustesUmeyama2d(transf, ptsA, ptsB);

    
    std::cout << "transfOut" << std::endl
              << transf.matrix() << std::endl;


    return 0;
}