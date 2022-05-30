/**
 * ARS - Angular Radon Spectrum 
 * Copyright (C) 2017-2020 Dario Lodi Rizzini
 *               2021- Dario Lodi Rizzini, Ernesto Fontana
 *
 * ARS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * ARS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with ARS.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <ars/ars2d.h>
#include <ars/ParamMap.h>
#include <ars/Profiler.h>
#include <ars/TranslationRefiner.h>

double mod180(double angle) {
    return (angle - floor(angle / M_PI) * M_PI);
}

//essential version of mpeg7 images dataset files reader/writer class
namespace ArsImgTests {

    class PointReaderWriter {
    public:
        ars::VectorVector2 points_;
        ars::Vector2 min_;
        ars::Vector2 max_;

        // Parameters read/written to files 
        double transl_x;
        double transl_y;
        double rot_theta;

        int num_in; //mpeg7 dataset related
        double noise_sigma;
        int num_occl;
        int num_rand;

        PointReaderWriter() {
        }

        PointReaderWriter(const PointReaderWriter& prw)
        : points_(prw.points_.begin(), prw.points_.end()), min_(prw.min_), max_(prw.max_) {
            transl_x = prw.transl_x;
            transl_y = prw.transl_y;
            rot_theta = prw.rot_theta;
            noise_sigma = prw.noise_sigma;
            num_in = prw.num_in;
            num_occl = prw.num_occl;
            num_rand = prw.num_rand;

            //            initFlann();
        }

        PointReaderWriter(const ars::VectorVector2& points)
        : points_(points.begin(), points.end()), min_(), max_() {
            transl_x = 0.0;
            transl_y = 0.0;
            rot_theta = 0.0;
            noise_sigma = 0.0;
            num_in = points.size();
            num_occl = points.size();
            num_rand = points.size();


        }

        virtual ~PointReaderWriter() {
        }

        int load(std::string filename) {
            ars::ParamMap params;
            //  char buffer[1000];
            std::string line, comment;
            ars::Vector2 p;
            size_t pos;
            int count;

            std::ifstream file(filename);
            if (!file) {
                std::cerr << "Cannot open file \"" << filename << "\"" << std::endl;
                return 0;
            }

            points_.clear();
            count = 0;
            while (!file.eof()) {
                std::getline(file, line);
                // Remove comments starting with '#'
                comment = "";
                pos = line.find_first_of('#');
                if (pos != std::string::npos) {
                    comment = line.substr(pos + 1, line.size());
                    line = line.substr(0, pos);
                }
                // Parse comment, if there is information
                std::stringstream sscomment(comment);
                params.read(sscomment);
                // Parse the line (after comment removal
                std::stringstream ssline(line);
                if (ssline >> p(0) >> p(1)) {
                    //      std::cout << "point [" << count << "]: " << p.transpose() << std::endl;
                    points_.push_back(p);
                    count++;
                }
            }
            file.close();

            // Computes coordinate bounds
            //        min_ << 1e+6, 1e+6;
            min_(0) = 1e+6;
            min_(1) = 1e+6;
            //        max_ << -1e+6, -1e+6;
            max_(0) = -1e+6;
            max_(1) = -1e+6;

            for (auto& p : points_) {
                if (p(0) < min_(0)) min_(0) = p(0);
                if (p(0) > max_(0)) max_(0) = p(0);
                if (p(1) < min_(1)) min_(1) = p(1);
                if (p(1) > max_(1)) max_(1) = p(1);
            }
            //  std::cout << "Read " << points_.size() << " points, min_ " << min_.transpose() << ", max_ " << max_.transpose() << std::endl;

            // Extracts params from ParamMap
            params.getParam<double>("transl_x", transl_x, 0.0);
            params.getParam<double>("transl_y", transl_y, 0.0);
            params.getParam<double>("rot_theta", rot_theta, 0.0);
            params.getParam<double>("noise_sigma", noise_sigma, noise_sigma);
            params.getParam<int>("num_in", num_in, count);
            params.getParam<int>("num_occl", num_occl, count);
            params.getParam<int>("num_rand", num_rand, count);
            //  std::cout << "\nParameters:\n";
            //  params.write(std::cout);

            //        // num_in is 0 at the beginning of the process
            //        if (num_in == 0) {
            //            num_in = count;
            //            std::cout << "**** EMPTY POINT SET ****" << std::endl;
            //        } else {
            //            initFlann();
            //        }

            return count;
        }

        ars::Vector2 computeCentroid() {
            ars::Vector2 mean = ars::Vector2::Zero();


            int count = 0;
            for (auto& p : points_) {
                mean += p;

                count++;
            }
            if (count > 0) {
                mean = mean / count;

            }
            return mean;
        }

        const ars::VectorVector2& points() const {
            return points_;
        }

        static Eigen::Affine2d coordToTransform(double x, double y, double theta) {
            Eigen::Affine2d transform = Eigen::Affine2d::Identity();


            transform.prerotate(theta);
            //        Eigen::preRotateAff2(transform, theta);
            transform.pretranslate(ars::Vector2(x, y));
            //        Eigen::preTranslateAff2(transform, x, y);

            return transform;
        }

        void applyTransform(const Eigen::Affine2d& transf) {
            min_ << 1e+6, 1e+6;
            max_ << -1e+6, -1e+6;

            for (auto& p : points_) {
                p = transf * p;
                //                ars::preTransfVec2(p, transf);

                if (p.x() < min_.x()) min_.x() = p.x();
                if (p.x() > max_.x()) max_.x() = p.x();
                if (p.y() < min_.y()) min_.y() = p.y();
                if (p.y() > max_.y()) max_.y() = p.y();
            }
            updateTransformInfo(transf);
        }

        void updateTransformInfo(const Eigen::Affine2d& transform) {
            Eigen::Affine2d prevTorig = coordToTransform(transl_x, transl_y, rot_theta);
            Eigen::Affine2d currTorig = transform * prevTorig;
            transformToCoodinates(currTorig, transl_x, transl_y, rot_theta);
        }

        void transformToCoodinates(const Eigen::Affine2d& transform, double& x, double& y, double& theta) {
            x = transform.matrix()(0, 2);
            //            x = transform.data_[0 * ars::Three + 2];
            y = transform.matrix()(1, 2);
            //            y = transform.data_[1 * ars::Three + 2];

            theta = atan2(transform.matrix()(1, 0), transform.matrix()(0, 0));
            //            theta = atan2(transform.data_[1 * ars::Three + 0], transform.data_[0 * ars::Three + 0]);

        }

        double getTranslX() const {
            return transl_x;
        }

        double getTranslY() const {
            return transl_y;
        }

        double getRotTheta() const {
            return rot_theta;
        }

        Eigen::Affine2d getTransform() const {
            Eigen::Translation2d t(transl_x, transl_y);
            Eigen::Affine3d rotBig = Eigen::Affine3d(Eigen::AngleAxisd(rot_theta, Eigen::Vector3d(0, 0, 1)));
            Eigen::Affine2d affMat = t * rotBig.linear().topLeftCorner<2, 2>();

            //        ARS_VAR4(affMat.matrix(), rot_theta, transl_x, transl_y);


            return affMat;
        }

    }; //end of class PointReaderWriter

} //end of namespace ArsImgTests

int main(int argc, char **argv) {
    ars::AngularRadonSpectrum2d arsSrc;
    ars::AngularRadonSpectrum2d arsDst;
    ArsImgTests::PointReaderWriter pointsSrc;
    ArsImgTests::PointReaderWriter pointsDst;


    ars::ParamMap params;
    std::string filenameCfg;
    std::string filenameSrc;
    std::string filenameDst;
    std::string filenameRot;
    std::string filenameArsSrc;
    std::string filenameArsDst;
    std::string filenameArsRot;
    std::string filenameArsCor;
    std::string filenameCovSrc;
    std::string filenameCovDst;
    int arsOrder;
    double arsSigma, arsThetaToll;
    double rotTrue, rotArs;
    //The variables below are for I/O related functionalities (plot, etc.) that are highly Eigen-based and are present in the CPU-only ArsImgTests...
    //Maybe implement them later
    //    double sampleRes, sampleAng; 
    //    int sampleNum;
    //    bool saveOn;
    //    bool saveCov;


    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, "");
    std::cout << "config filename: " << filenameCfg << std::endl;
    if (filenameCfg != "") {
        params.read(filenameCfg);
    }

    params.read(argc, argv);
    params.getParam<std::string>("src", filenameSrc, "/home/rimlab/Downloads/mpeg7_point_tests/noise000_occl00_rand000/apple-1_xp0686_yp0967_t059_sigma0001_occl000.txt");
    params.getParam<std::string>("dst", filenameDst, "/home/rimlab/Downloads/mpeg7_point_tests/noise000_occl00_rand000/apple-1_xp0749_yn0521_t090_sigma0001_occl000.txt");
    params.getParam<int>("arsOrder", arsOrder, 20);
    params.getParam<double>("arsSigma", arsSigma, 1.0);
    params.getParam<double>("arsTollDeg", arsThetaToll, 1.0);
    arsThetaToll *= M_PI / 180.0;
    //    params.getParam<double>("sampleResDeg", sampleRes, 0.5);
    //    sampleRes *= M_PI / 180.0;
    //    params.getParam<bool>("saveOn", saveOn, false);
    //    params.getParam<bool>("saveCov", saveCov, false);

    std::cout << "\nParameter values:\n";
    params.write(std::cout);
    std::cout << std::endl;


    // Loads files and computes the rotation
    std::cout << "\n*****\nLoading file \"" << filenameSrc << "\"" << std::endl;
    pointsSrc.load(filenameSrc);
    std::cout << "\n*****\nLoading file \"" << filenameDst << "\"" << std::endl;
    pointsDst.load(filenameDst);
    std::cout << "  points src " << pointsSrc.points().size() << ", points dst " << pointsDst.points().size() << std::endl;

    int numPts = std::min<int>(pointsSrc.points().size(), pointsDst.points().size()); //the two should normally be equals
    //    int numPtsAfterPadding = numPts;

    //ARS parameters setting
    arsSrc.setARSFOrder(arsOrder);
    arsDst.setARSFOrder(arsOrder);
    ars::ArsKernelIsotropic2d::ComputeMode pnebiMode = ars::ArsKernelIsotropic2d::ComputeMode::PNEBI_DOWNWARD;
    arsSrc.setComputeMode(pnebiMode);
    arsDst.setComputeMode(pnebiMode);


    std::cout << "\n------\n" << std::endl;





    //ARS SRC 

    //END OF ARS SRC



    std::cout << "\n------\n" << std::endl; //"pause" between ars src and ars dst



    //ARS DST -> preparation for kernel calls and kernel calls


    //END OF ARS DST





    //Computation final computations (correlation, ...) on CPU
    //    std::cout << "\nARS Coefficients:\n";
    //    std::cout << "Coefficients: Src, Dst, Cor" << std::endl;
    //
    //    double thetaMax, corrMax, fourierTol;
    //    fourierTol = 1.0; // TODO: check for a proper tolerance
    //
    //    std::vector<double> coeffsCor;
    //    {
    //        ars::ScopedTimer("ars.correlation()");
    //        std::vector<double> tmpSrc;
    //        tmpSrc.assign(coeffsArsSrc, coeffsArsSrc + coeffsMatNumColsPadded);
    //        std::vector<double> tmpDst;
    //        tmpDst.assign(coeffsArsDst, coeffsArsDst + coeffsMatNumColsPadded);
    //        ars::computeFourierCorr(tmpSrc, tmpDst, coeffsCor);
    //        ars::findGlobalMaxBBFourier(coeffsCor, 0.0, M_PI, arsThetaToll, fourierTol, thetaMax, corrMax);
    //        rotArs = thetaMax;
    //    }
    //
    //
    //
    //    arsSrc.setCoefficients(coeffsArsSrc, 2 *);
    //    //    for (int i = 0; i < coeffsVectorMaxSz; i++) {
    //    //        std::cout << "arsSrc - coeff_d[" << i << "] " << d_coeffsMat1[i] << std::endl;
    //    //    }
    //    arsDst.setCoefficients(coeffsArsDst, coeffsMatNumCols);
    //    for (int i = 0; i < arsSrc.coefficients().size() && i < arsDst.coefficients().size(); ++i) {
    //        std::cout << "\t" << i << " \t" << arsSrc.coefficients().at(i) << " \t" << arsDst.coefficients().at(i) << " \t" << coeffsCor[i] << std::endl;
    //    }
    //    std::cout << std::endl;



    // Computes the rotated points, centroid, affine transf matrix between src and dst
    ArsImgTests::PointReaderWriter pointsRot(pointsSrc.points());
    ars::Vector2 centroidSrc = pointsSrc.computeCentroid();
    ars::Vector2 centroidDst = pointsDst.computeCentroid();
    Eigen::Affine2d rotSrcDst = ArsImgTests::PointReaderWriter::coordToTransform(0.0, 0.0, rotArs);
    Eigen::Translation2d translSrcDst(centroidDst - rotSrcDst * centroidSrc);
    Eigen::Affine2d transfSrcDst = rotSrcDst * translSrcDst;

    //    std::cout << "centroidSrc " << centroidSrc.transpose() << "\n"
    //            << "rotSrcDst\n" << rotSrcDst.matrix() << "\n"
    //            << "translation: [" << translSrcDst.transpose() << "] rotation[deg] " << (180.0 / M_PI * rotArs) << "\n";
    std::cout << "centroidSrc " << centroidSrc(0) << " \t" << centroidSrc(1) << "\n"
            << "centroidDst " << centroidDst(0) << " \t" << centroidDst(1) << "\n"
            << "rotSrcDst\n" << rotSrcDst.matrix() << "\n"
            << "translation: [" << translSrcDst.x() << " \t" << translSrcDst.y() << "] rotation[deg] " << (180.0 / M_PI * rotArs) << "\n";

    pointsRot.applyTransform(transfSrcDst);



    rotTrue = pointsDst.getRotTheta() - pointsSrc.getRotTheta();
    std::cout << "\n***\npointsDst.getrotTheta() [deg]" << (180 / M_PI * pointsDst.getRotTheta())
            << ", pointsSrc.getrotTheta() [deg] " << (180.0 / M_PI * pointsSrc.getRotTheta()) << "\n";
    std::cout << "rotTrue[deg] \t" << (180.0 / M_PI * rotTrue) << " \t" << (180.0 / M_PI * mod180(rotTrue)) << std::endl;
    std::cout << "rotArs[deg] \t" << (180.0 / M_PI * rotArs) << " \t" << (180.0 / M_PI * mod180(rotArs)) << std::endl;




    return 0;
}

