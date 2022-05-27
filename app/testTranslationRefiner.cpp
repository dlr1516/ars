#include <iostream>

#include <ars/TranslationRefiner.h>

int main(int argc, char **argv) {

    std::cout << "This program executes a simple version of the standard registration algorithm Iterative Closest Point (ICP)" << std::endl;

    // define a 3 dim problem with 10000 model points
    // and 10000 template points:
    int32_t dim = 2;
    int32_t num = 100;


    ars::VectorVector2 modelPts, templatePts;

    // set model and template points
    std::cout << std::endl << "Creating model with 10000 points ..." << std::endl;
    std::cout << "Creating template by shifting model by (2,1) ..." << std::endl;
    int32_t k = 0;
    for (double x = -2; x < 2; x += 0.04) {
        double y = 5 * x * exp(-x * x);
        ars::Vector2 ptM(x, y);
        modelPts.push_back(ptM);
        ars::Vector2 ptT(x - 2, y - 1);
        templatePts.push_back(ptT);

        k++;
        //        std::std::cout << "k " << k << std::std::endl;
    }

    // version with Eigen
    Eigen::Affine2d transf;
    transf.setIdentity();

    ars::TranslationRefiner translRefiner(modelPts, templatePts, transf);
    translRefiner.computeProcrustes(transf);
    std::cout << "transf:" << std::endl << transf.matrix() << std::endl;


    return 0;
}
