#ifndef POINT_READER_WRITER_H
#define POINT_READER_WRITER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <random>


#include <flann/flann.hpp>
#include <flann/flann.h>
//#include <pcl/kdtree/kdtree_flann.h>
//
//#include <pcl/io/pcd_io.h>
//#include <pcl/point_cloud.h> 
//#include <pcl/point_types.h>

#include "definitions.h"
#include "utils.h"

#include <rofl/common/param_map.h>

#define POINT_READER_WRITER_GET(TYPE,X) TYPE get##X() const { return (X); }

namespace ArsImgTests {

    std::pair<double, double> addInverval(const std::pair<double, double>& interv1, const std::pair<double, double>& interv2);

    /** Class to read, write and store point sets. 
     */
    class PointReaderWriter {
    public:


        /** Constructor.
         */
        PointReaderWriter();

        /** Contructor with the possibility to read parameters from 
         */
        PointReaderWriter(std::istream& in);

        /** Constructor.
         */
        PointReaderWriter(const PointReaderWriter& prw);

        /** Constructor.
         */
        PointReaderWriter(const cuars::VecVec2d& points);
        /**
         */
        PointReaderWriter(const cuars::VecVec2d& points, const cuars::Affine2d& transf);

        /**
         * Destructor.
         */
        ~PointReaderWriter();

        /**
         */
        void swap(PointReaderWriter& prw);

        /** Loads points (and potential parameters) from file.
         * If some parameters are defined in the comments of the point file, 
         * the current parameters are overwritten. 
         */
        int load(std::string filename);

        /** Loads points (and potential parameters) from pcd file.
         */
        //        void loadCloudPoints(std::string cloudFilename);

        /** Saves points (and potential parameters) from file.
         */
        void save(std::string filename);

        /**
         */
        const cuars::VecVec2d& points() const;

        void reset();

        void insertPoints(const cuars::VecVec2d& points);

        double xmin() const;

        double xmax() const;

        double ymin() const;

        double ymax() const;

        std::pair<double, double> xinterv() const;

        std::pair<double, double> yinterv() const;

        //        POINT_READER_WRITER_GET(double, transl_x)
        //        POINT_READER_WRITER_GET(double, transl_y)
        //        POINT_READER_WRITER_GET(double, rot_theta)
        //        POINT_READER_WRITER_GET(double, noise_sigma)
        //        POINT_READER_WRITER_GET(int, num_in)
        //        POINT_READER_WRITER_GET(int, num_occl)
        //        POINT_READER_WRITER_GET(int, num_rand)

        double getTranslX() const;

        double getTranslY() const;

        double getRotTheta() const;

        double getNoiseSigma() const;

        int getNumIn() const;

        int getNumOccl() const;

        int getNumRand() const;


        // --------------------------------------------------------
        // POINT TRANSFORMATION METHODS
        // --------------------------------------------------------

        void findKNearest(const cuars::Vec2d& query, int k, std::vector<int>& indices, std::vector<double>& distances) const;

        void computeContour(cuars::VecVec2d& contour) const;

        // --------------------------------------------------------
        // POINT TRANSFORMATION METHODS
        // --------------------------------------------------------

        /** Creates random transformation and applies it to points.
         * Parameters of random transformation are computed.
         */
        void applyRandTransform(double translationMax, double rotMin, double rotMax);

        /** Applies a transformation to all the points.
         * It composes transformation with previous transformation. 
         */
        void applyTransform(const cuars::Affine2d& transf);

        /** Applies a transformation to all the points.
         * It composes transformation with previous transformation. 
         */
        void applyTransform(double x, double y, double theta);

        /** Applies occlusion to points for a given occlusion percentage (this percentage 
         * is only a desired target). 
         */
        void applyOcclusion(double occlPerc = 0.15);

        /** Adds random (gaussian) noise to the points.
         */
        void addRandNoise(double noiseSigma);

        /** Adds new random points to the set (percentage of the current set).
         */
        void addRandPoints(double perc = 0.10, int maxnum = -1);

        /** Removes bias, i.e. it subtracts the mean point to all points coordinates.
         */
        void removeBias();

        /** Computes the centroid.
         */
        cuars::Vec2d computeCentroid() const;

        /** Computes the covariance.
         */
        cuars::Mat2d computeCovariance() const;

        // --------------------------------------------------------
        // STATIC METHODS
        // --------------------------------------------------------

        /** Coordinates [x,y,theta] to transform matrix
         */
        static cuars::Affine2d coordToTransform(double x, double y, double theta);

        /** Transform matrix to coordinates.
         */
        void transformToCoodinates(const cuars::Affine2d& transform, double& x, double& y, double& theta);

    private:
        cuars::VecVec2d points_;
        cuars::Vec2d min_;
        cuars::Vec2d max_;
        std::random_device randDev_;
        std::mt19937 randGen_;

        // Parameters read/written to files 
        double transl_x;
        double transl_y;
        double rot_theta;
        double noise_sigma;
        int num_in;
        int num_occl;
        int num_rand;

        // Parameters for search
        flann::Index<flann::L2<double> >* flannIndexer_;
        flann::Matrix<double> flannPoints_;
        std::vector<double> pointsVec_;

        /** Updates the data about the random transformation by composing previous transform 
         * with the given one. 
         */
        void updateTransformInfo(const cuars::Affine2d& transform);

        /**
         * Creates the FLANN data structure. 
         */
        void initFlann();

        /**
         * Deallocates FLANN data structure. 
         */
        void clearFlann();

        /**
         */
        static bool ccw(const cuars::Vec2d& a, const cuars::Vec2d& b, const cuars::Vec2d& c) {
            cuars::Vec2d ab, ac;
            cuars::vec2diff(ab, b, a);
            cuars::vec2diff(ac, c, a);

            return (ab.x * ac.y - ab.y * ac.x);
        }
    };

} // end of namespace 

#endif