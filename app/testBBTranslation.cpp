#include <iostream>

#include <ars/BBTranslation.h>

void findBoundingBox(const ars::VectorVector2& pts,
                     ars::Vector2& ptMin,
                     ars::Vector2& ptMax);

int main(int argc, char** argv) {
    ars::VectorVector2 ptsA, ptsB;
    ars::Vector2 minA, maxA, minB, maxB;
    int dim = 3;
    Eigen::Affine2d transf;
    ars::BBTranslation translEstim;

    Eigen::Vector2d aA(Eigen::Vector2d(1, 1));
    Eigen::Vector2d bA(Eigen::Vector2d(3, 4));
    Eigen::Vector2d cA(Eigen::Vector2d(10, 8));
    Eigen::Vector2d dA(Eigen::Vector2d(4, 7));
    Eigen::Vector2d eA(Eigen::Vector2d(7, 2));
    Eigen::Vector2d fA(Eigen::Vector2d(-14, -5));
    Eigen::Vector2d zeroB(Eigen::Vector2d(0, 0));
    ptsA.push_back(aA);
    ptsA.push_back(bA);
    ptsA.push_back(cA);
    ptsA.push_back(dA);
    ptsA.push_back(eA);
    ptsA.push_back(fA);
    Eigen::Affine2d transfTrue;
    // ars::createRotationMatrix(transfTrue, 0, 0, 30 * M_PI / 180.0);
    transfTrue.translation() = Eigen::Vector2d(-1, -1);
    std::cout << "Applying transformation transfTrue:" << std::endl
              << transfTrue.matrix() << std::endl;
    for (int i = 0; i < 6; ++i)  // 6 = number of points
    {
        Eigen::Vector2d ptsAiTransf =
            ptsA[i] + transfTrue.translation().matrix();
        ptsB.push_back(Eigen::Vector2d(ptsAiTransf(0), ptsAiTransf(1)));
    }

    findBoundingBox(ptsA, minA, maxA);
    findBoundingBox(ptsB, minB, maxB);
    std::cout << "ptsA: size " << ptsA.size() << ", min [" << minA.transpose()
              << "]  max [" << maxA.transpose() << "]" << std::endl;
    std::cout << "ptsB: size " << ptsB.size() << ", min [" << minB.transpose()
              << "]  max [" << maxB.transpose() << "]" << std::endl;
    std::cout << "translation interval: min [" << (minB - maxA).transpose()
              << "]  max [" << (maxB - minA).transpose() << "]" << std::endl;

    translEstim.setTranslMinMax(minB - maxA, maxB - minA);
    translEstim.setPts(ptsA, ptsB);
    translEstim.compute();

    return 0;
}

void findBoundingBox(const ars::VectorVector2& pts,
                     ars::Vector2& ptMin,
                     ars::Vector2& ptMax) {
    for (int i = 0; i < pts.size(); ++i) {
        for (int d = 0; d < 2; ++d) {
            if (i == 0 || pts[i](d) < ptMin(d)) {
                ptMin(d) = pts[i](d);
            }
            if (i == 0 || pts[i](d) > ptMax(d)) {
                ptMax(d) = pts[i](d);
            }
        }
    }
}