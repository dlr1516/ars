#include <algorithm>


#include "ars/mpeg7RW.h"


namespace ArsImgTests {

    std::pair<double, double> addInverval(const std::pair<double, double>& interv1, const std::pair<double, double>& interv2) {
        std::pair<double, double> intervAdd = std::make_pair(std::min(interv1.first, interv2.first), std::max(interv2.second, interv2.second));
        return intervAdd;
    }

    PointReaderWriter::PointReaderWriter()
    : points_(), min_(), max_(), randDev_(), randGen_(randDev_()), //
    flannIndexer_(nullptr), flannPoints_(), pointsVec_() {
        transl_x = 0.0;
        transl_y = 0.0;
        rot_theta = 0.0;
        noise_sigma = 1.0;
        num_in = 0;
        num_occl = 0;
        num_rand = 0;
    }

    PointReaderWriter::PointReaderWriter(std::istream& in)
    : points_(), min_(), max_(), randDev_(), randGen_(randDev_()),
    flannIndexer_(nullptr), flannPoints_(), pointsVec_() {
        rofl::ParamMap params;
        params.read(in);
        params.getParam<double>("transl_x", transl_x, 0.0);
        params.getParam<double>("transl_y", transl_y, 0.0);
        params.getParam<double>("rot_theta", rot_theta, 0.0);
        params.getParam<double>("noise_sigma", noise_sigma, 1.0);
        params.getParam<int>("num_in", num_in, 0);
        params.getParam<int>("num_occl", num_occl, 0);
        params.getParam<int>("num_rand", num_rand, 0);
    }

    PointReaderWriter::PointReaderWriter(const PointReaderWriter& prw)
    : points_(prw.points_.begin(), prw.points_.end()), min_(prw.min_), max_(prw.max_), randDev_(), randGen_(randDev_()),
    flannIndexer_(nullptr), flannPoints_(), pointsVec_() {
        transl_x = prw.transl_x;
        transl_y = prw.transl_y;
        rot_theta = prw.rot_theta;
        noise_sigma = prw.noise_sigma;
        num_in = prw.num_in;
        num_occl = prw.num_occl;
        num_rand = prw.num_rand;

        initFlann();
    }

    PointReaderWriter::PointReaderWriter(const cuars::VecVec2d& points)
    : points_(points.begin(), points.end()), min_(), max_(), randDev_(), randGen_(randDev_()),
    flannIndexer_(nullptr), flannPoints_(), pointsVec_() {
        transl_x = 0.0;
        transl_y = 0.0;
        rot_theta = 0.0;
        noise_sigma = 0.0;
        num_in = points.size();
        num_occl = points.size();
        num_rand = points.size();

        initFlann();
    }

    PointReaderWriter::PointReaderWriter(const cuars::VecVec2d& points, const cuars::Affine2d& transf)
    : points_(points.begin(), points.end()), min_(), max_(), randDev_(), randGen_(randDev_()),
    flannIndexer_(nullptr), flannPoints_(), pointsVec_() {
        transformToCoodinates(transf, transl_x, transl_y, rot_theta);
        noise_sigma = 0.0;
        num_in = points.size();
        num_occl = points.size();
        num_rand = points.size();

        //  std::cout << "Init with transform: transl [" << transl_x << "," << transl_y << "], rot " << rot_theta << "\n" << transf.matrix() << std::endl;

        for (auto& p : points_) {
            //    std::cout << "  " << p.transpose() << " --> ";
            p = transf * p;
            //    std::cout << p.transpose() << std::endl;
        }

        initFlann();
    }

    //PointReaderWriter::PointReaderWriter(const std::vector<cuars::Vec2d>& points,const cuars::Affine2d& transf)
    //  : points_(), min_(), max_(), randDev_(), randGen_(randDev_()) 
    //{
    //  transformToCoodinates(transf,transl_x,transl_y,rot_theta);
    //  noise_sigma = 0.0;
    //  num_in = points.size();
    //  num_occl = points.size();
    //  num_rand = points.size();

    //  std::cout << "Init with transform: transl [" << transl_x << "," << transl_y << "], rot " << rot_theta << "\n" << transf.matrix() << std::endl;

    //  for (auto& p : points) {
    //    points_.push_back(transf * p);
    //  }
    //}

    PointReaderWriter::~PointReaderWriter() {
        clearFlann();
    }

    void PointReaderWriter::swap(PointReaderWriter& prw) {
        std::swap(points_, prw.points_);
        std::swap(min_, prw.min_);
        std::swap(max_, prw.max_);
        std::swap(transl_x, prw.transl_x);
        std::swap(transl_y, prw.transl_y);
        std::swap(rot_theta, prw.rot_theta);
        std::swap(num_in, prw.num_in);
        std::swap(num_occl, prw.num_occl);
        std::swap(num_rand, prw.num_rand);
    }

    int PointReaderWriter::load(std::string filename) {
        rofl::ParamMap params;
        //  char buffer[1000];
        std::string line, comment;
        cuars::Vec2d p;
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
            // Parse comment, if there are information
            std::stringstream sscomment(comment);
            params.read(sscomment);
            // Parse the line (after comment removal
            std::stringstream ssline(line);
            if (ssline >> p.x >> p.y) {
                //      std::cout << "point [" << count << "]: " << p.transpose() << std::endl;
                points_.push_back(p);
                count++;
            }
        }
        file.close();

        // Computes coordinate bounds
        //        min_ << 1e+6, 1e+6;
        min_.x = 1e+6;
        min_.y = 1e+6;
        //        max_ << -1e+6, -1e+6;
        max_.x = -1e+6;
        max_.y = -1e+6;

        for (auto& p : points_) {
            if (p.x < min_.x) min_.x = p.x;
            if (p.x > max_.x) max_.x = p.x;
            if (p.y < min_.y) min_.y = p.y;
            if (p.y > max_.y) max_.y = p.y;
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

        // num_in is 0 at the beginning of the process
        if (num_in == 0) {
            num_in = count;
            std::cout << "**** EMPTY POINT SET ****" << std::endl;
        } else {
            initFlann();
        }

        return count;
    }

    //    void PointReaderWriter::loadCloudPoints(std::string cloudFilename) {
    //        // Loads the input point cloud
    //        points_.clear();
    //        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    //
    //        std::cout << "Loading point cloud from \"" << cloudFilename << "\"" << std::endl;
    //        if (pcl::io::loadPCDFile(cloudFilename, *cloud) < 0) {
    //            std::cerr << "Cannot load point cloud from \"" << cloudFilename << "\"" << std::endl;
    //            return;
    //        }
    //
    //        size_t cloudSz = cloud->size();
    //        for (size_t i = 0; i < cloudSz; ++i) {
    //            cuars::Vec2d ptEig;
    //            ptEig(0) = cloud->points.at(i).x;
    //            ptEig(1) = cloud->points.at(i).y;
    //
    //            points_.push_back(ptEig);
    //
    //        }
    //        //computing rotation on z-axis and assigning it to rot_theta
    //        Eigen::Rotation2D<float> rotMat(cloud->sensor_orientation_.toRotationMatrix().block<2, 2>(0, 0));
    //        rot_theta = rotMat.angle(); //using it in radians during the algorithm
    //    }

    void PointReaderWriter::save(std::string filename) {
        std::ofstream file(filename);
        if (!file) {
            std::cerr << "Cannot open file \"" << filename << "\"" << std::endl;
            return;
        }
        file << "# transl_x " << transl_x << "\n"
                << "# transl_y " << transl_y << "\n"
                << "# rot_theta " << rot_theta << "\n"
                << "# noise_sigma " << noise_sigma << "\n"
                << "# num_in " << num_in << "\n"
                << "# num_occl " << num_occl << "\n"
                << "# num_rand " << num_rand << "\n";
        for (auto& p : points_) {
            file << p.x << " " << p.y << "\n";
        }
        file.close();
    }

    /**
     */
    const cuars::VecVec2d& PointReaderWriter::points() const {
        return points_;
    }

    void PointReaderWriter::reset() {
        if (points_.size() == 0)
            return;
        points_.clear();
        clearFlann();
    }

    void PointReaderWriter::insertPoints(const cuars::VecVec2d& points) {
        points_ = points;
        initFlann();
    }

    double PointReaderWriter::xmin() const {
        return min_.x;
    }

    double PointReaderWriter::xmax() const {
        return max_.x;
    }

    double PointReaderWriter::ymin() const {
        return min_.y;
    }

    double PointReaderWriter::ymax() const {
        return max_.y;
    }

    std::pair<double, double> PointReaderWriter::xinterv() const {
        return std::make_pair(min_.x, max_.x);
    }

    std::pair<double, double> PointReaderWriter::yinterv() const {
        return std::make_pair(min_.y, max_.y);
    }

    //        POINT_READER_WRITER_GET(double, transl_x)
    //        POINT_READER_WRITER_GET(double, transl_y)
    //        POINT_READER_WRITER_GET(double, rot_theta)
    //        POINT_READER_WRITER_GET(double, noise_sigma)
    //        POINT_READER_WRITER_GET(int, num_in)
    //        POINT_READER_WRITER_GET(int, num_occl)
    //        POINT_READER_WRITER_GET(int, num_rand)

    double PointReaderWriter::getTranslX() const {
        return transl_x;
    }

    double PointReaderWriter::getTranslY() const {
        return transl_y;
    }

    double PointReaderWriter::getRotTheta() const {
        return rot_theta;
    }

    double PointReaderWriter::getNoiseSigma() const {
        return noise_sigma;
    }

    int PointReaderWriter::getNumIn() const {
        return num_in;
    }

    int PointReaderWriter::getNumOccl() const {
        return num_occl;
    }

    int PointReaderWriter::getNumRand() const {
        return num_rand;
    }

    void PointReaderWriter::findKNearest(const cuars::Vec2d& query, int k, std::vector<int>& indices, std::vector<double>& distances) const {
        double queryVec[2];

        // Initialization of query FLANN matrix
        queryVec[0] = query.x;
        queryVec[1] = query.y;
        flann::Matrix<double> flannQuery = flann::Matrix<double>(queryVec, 1, 2);
        // Finds the k-nearest of current point 
        indices.resize(k, -1);
        flann::Matrix<int> flannIndices(&indices[0], 1, (int) k);
        distances.resize(k);
        flann::Matrix<double> flannDistances(&distances[0], 1, (int) k);

        // Calls FLANN search method
        flann::SearchParams param;
        param.eps = 1.0;
        param.sorted = true;
        param.checks = 128;
        // flannIndexer_ may be null with empty point clouds
        if (flannIndexer_ != nullptr) {
            flannIndexer_->knnSearch(flannQuery, flannIndices, flannDistances, k, param);
        }
    }

    void PointReaderWriter::computeContour(cuars::VecVec2d& contour) const {
        std::vector<bool> visited(points_.size(), false);
        std::vector<int> indices;
        std::vector<double> distances;
        cuars::Vec2d center = computeCentroid();
        int k = 50;
        int idx, inext;

        // Computes the contour
        contour.clear();
        idx = 0;
        while (idx >= 0) {
            if (0 <= idx && idx < points_.size() && visited[idx] == false) {
                visited[idx] = true;
                contour.push_back(points_[idx]);
                findKNearest(points_[idx], k, indices, distances);
                inext = -1;
                for (auto& i : indices) {
                    if (0 <= i && i < points_.size() && visited[i] == false && ccw(center, points_[idx], points_[i])) {
                        inext = i;
                        break;
                    }
                }
                idx = inext;
            }
        }
    }

    void PointReaderWriter::applyRandTransform(double translationMax, double rotMin, double rotMax) {
        std::uniform_real_distribution<> randX(-translationMax, translationMax);
        std::uniform_real_distribution<> randY(-translationMax, translationMax);
        std::uniform_real_distribution<> randTheta(rotMin, rotMax);
        applyTransform(randX(randGen_), randY(randGen_), randTheta(randGen_));
    }

    void PointReaderWriter::applyTransform(const cuars::Affine2d& transf) {
        //        min_ << 1e+6, 1e+6;
        min_.x = 1e+6;
        min_.y = 1e+6;
        //        max_ << -1e+6, -1e+6;
        max_.x = -1e+6;
        max_.y = -1e+6;

        for (auto& p : points_) {
            p = transf * p;
            if (p.x < min_.x) min_.x = p.x;
            if (p.x > max_.x) max_.x = p.x;
            if (p.y < min_.y) min_.y = p.y;
            if (p.y > max_.y) max_.y = p.y;
        }
        updateTransformInfo(transf);
    }

    void PointReaderWriter::applyTransform(double x, double y, double theta) {
        std::cout << __FILE__ << "," << __LINE__ << ": transform x " << x << ", y " << y << ", theta[deg] " << (180.0 / M_PI * theta) << std::endl;
        cuars::Affine2d transf = coordToTransform(x, y, theta);
        applyTransform(transf);
    }

    void PointReaderWriter::applyOcclusion(double occlPerc) {
        double sx = max_.x - min_.x;
        double sy = max_.y - min_.y;
        double radiusMean = occlPerc * sqrt(sx * sy);
        double radiusSigma = 0.05 * radiusMean;
        //  std::uniform_real_distribution<> randX(min_.x,max_.x);
        //  std::uniform_real_distribution<> randY(min_.y,max_.y);
        std::uniform_int_distribution<> randIdx(0, points_.size() - 1);

        // Generate occlusion region as a random circle 
        //  cuars::Vec2d center(randX(randGen_),randY(randGen_));
        //  double radius = randRadius(randGen_);
        cuars::Vec2d center = points_[randIdx(randGen_)];
        double radius = radiusMean;
        cuars::VecVec2d pointsTmp;
        for (auto& p : points_) {
            if ((p - center).norm() > radius) {
                pointsTmp.push_back(p);
            }
        }
        std::swap(points_, pointsTmp);
        num_occl = points_.size();
        // We should update the bounds min_ and max_, but since the new bounds are 
        // inside this operation is skipped.
    }

    void PointReaderWriter::addRandNoise(double noiseSigma) {
        std::normal_distribution<> randNoise(0.0, noiseSigma);
        for (auto& p : points_) {
            p.x += randNoise(randGen_);
            p.y += randNoise(randGen_);
        }
        noise_sigma = noiseSigma; // parameter saved with the point file
    }

    void PointReaderWriter::addRandPoints(double perc, int maxnum) {
        double dx = max_.x - min_.x;
        double dy = max_.y - min_.y;
        //  std::uniform_real_distribution<> randX(min_.x-dx,max_.x+dx);
        //  std::uniform_real_distribution<> randY(min_.y-dy,max_.y+dy);
        //  int num = perc * points_.size();
        //  cuars::Vec2d p;
        //  for (int i = 0; i < num; ++i) {
        //    p << randX(randGen_), randY(randGen_);
        //    points_.push_back(p);
        //  }
        std::uniform_real_distribution<> randRadius(0.0, std::max(dx, dy));
        std::uniform_real_distribution<> randTheta(0.0, 2.0 * M_PI);
        int num = perc * points_.size();
        cuars::Vec2d p;
        double r, a;
        for (int i = 0; i < num; ++i) {
            r = randRadius(randGen_);
            a = randRadius(randGen_);
            //            p << r * cos(a), r * sin(a);
            p.x = r * cos(a);
            p.y = r * sin(a);

            //            p += 0.5 * (min_ + max_);
            p.x += 0.5 * (min_.x + max_.x);
            p.y += 0.5 * (min_.y + max_.y);

            points_.push_back(p);
        }

        // If there is a maximum number of points (maxnum >= 0), then the whole set is resampled.
        if (maxnum >= 0 && points_.size() >= maxnum) {
            std::uniform_int_distribution<> randIdx(0, points_.size() - 1);
            for (int i = 0; i < maxnum; ++i) {
                std::swap(points_[i], points_[ randIdx(randGen_) ]);
            }
            points_.erase(points_.begin() + maxnum, points_.end());
        }

        num_rand = points_.size();
    }

    void PointReaderWriter::removeBias() {
        //        cuars::Vec2d mean = cuars::Vec2d::Zero();
        cuars::Vec2d mean;
        mean.x = 0.0;
        mean.y = 0.0;
        int count = 0;
        for (auto& p : points_) {
            if (p.x < min_.x) min_.x = p.x;
            if (p.x > max_.x) max_.x = p.x;
            if (p.y < min_.y) min_.y = p.y;
            if (p.y > max_.y) max_.y = p.y;
            mean.x += p.x;
            mean.y += p.y;

            count++;
        }
        if (count > 0) {
            mean.x /= count;
            mean.y /= count;
        }
        for (auto& p : points_) {
            p.x -= mean.x;
            p.y -= mean.y;
        }
        //        min_ -= mean;
        min_.x -= mean.x;
        min_.y -= mean.y;
        //                max_ -= mean;
        max_.x -= mean.x;
        max_.y -= mean.y;

    }

    cuars::Vec2d PointReaderWriter::computeCentroid() const {
        cuars::Vec2d mean;
        mean.x = 0.0;
        mean.y = 0.0;

        int count = 0;
        for (auto& p : points_) {
            mean.x += p.x; //maybe use cuars::vec2sum function instead of going element by element?
            mean.y += p.y; //maybe use cuars::vec2sum function instead of going element by element?

            count++;
        }
        if (count > 0) {
            //            mean = mean / count;
            mean.x /= count;
            mean.y /= count;

        }
        return mean;
    }

    cuars::Mat2d PointReaderWriter::computeCovariance() const {
        cuars::Vec2d mean = computeCentroid();
        //        cuars::Mat2d cov = cuars::Mat2d::Zero();
        cuars::Mat2d cov;
        cov.w = 0.0;
        cov.x = 0.0;
        cov.y = 0.0;
        cov.z = 0.0;

        int count = 0;
        for (auto& p : points_) {
            cov += (p - mean) * (p - mean).transpose();
            count++;
        }
        if (count > 0) {
            //            cov = cov / count;
            cov.w /= count;
            cov.x /= count;
            cov.y /= count;
            cov.z /= count;



        }
        return cov;
    }


    // --------------------------------------------------------
    // STATIC METHODS
    // --------------------------------------------------------

    cuars::Affine2d PointReaderWriter::coordToTransform(double x, double y, double theta) {
        //        cuars::Affine2d transform = cuars::Affine2d::Identity();
        cuars::Affine2d transform(0.0, 0.0, 0.0);


        transform.prerotate(theta);
        transform.pretranslate(cuars::Vec2d(x, y));
        return transform;
    }

    void PointReaderWriter::transformToCoodinates(const cuars::Affine2d& transform, double& x, double& y, double& theta) {
        //        x = transform.matrix()(0, 2);
        x = transform.data_[0 * cuars::Three + 2];
        //        y = transform.matrix()(1, 2);
        y = transform.data_[1 * cuars::Three + 2];

        //        theta = atan2(transform.matrix()(1, 0), transform.matrix()(0, 0));
        theta = atan2(transform.data_[1 * cuars::Three + 0], transform.data_[0 * cuars::Three + 0]);

    }

    // --------------------------------------------------------
    // PRIVATE METHODS
    // --------------------------------------------------------

    void PointReaderWriter::updateTransformInfo(const cuars::Affine2d& transform) {
        cuars::Affine2d prevTorig = coordToTransform(transl_x, transl_y, rot_theta);
        cuars::Affine2d currTorig = transform * prevTorig;
        transformToCoodinates(currTorig, transl_x, transl_y, rot_theta);
    }

    void PointReaderWriter::initFlann() {
        if (points_.empty()) {
            flannIndexer_ = nullptr;
        }
        pointsVec_.resize(2 * points_.size());
        for (int i = 0; i < points_.size(); ++i) {
            pointsVec_[2 * i] = points_[i].x;
            pointsVec_[2 * i + 1] = points_[i].y;
        }
        flannPoints_ = flann::Matrix<double>(&pointsVec_[0], points_.size(), 2);
        flannIndexer_ = new flann::Index<flann::L2<double> >(flannPoints_, flann::KDTreeIndexParams(4));
        flannIndexer_->buildIndex();
    }

    void PointReaderWriter::clearFlann() {
        if (flannIndexer_ != 0) {
            delete flannIndexer_;
        }
    }

} // end of namespace 
