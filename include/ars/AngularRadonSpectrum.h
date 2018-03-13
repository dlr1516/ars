#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <ars/functions.h>

namespace ars
{

// --------------------------------------------------------
// ARSF FUNCTIONS
// --------------------------------------------------------

/** Computes the contribution of term (lambda,phi) to the Angular Radon Spectrum Fourier (ARSF)
 * coefficients of ARS.
 */
void updateARSFCoeffDirect(double lambda,double phi,int n,Eigen::VectorXd& coeffs);

/** Computes the contribution of term (lambda,phi) to the Angular Radon Spectrum Fourier (ARSF) 
 * coefficients computed with upward recursion of modified Bessel functions.
 * WARNING: Upward recursion is numerically UNSTABLE!
 */
void updateARSFCoeffRecursUp(double lambda,double phi,int n,Eigen::VectorXd& coeffs);

/** Computes the contribution of term (lambda,phi) to the Angular Radon Spectrum Fourier (ARSF) 
 * coefficients computed with downward recursion of modified Bessel functions.
 */
void updateARSFCoeffRecursDown(double lambda,double phi,int n,Eigen::VectorXd& coeffs);

/** Computes the contribution of term (lambda,phi) to the Angular Radon Spectrum Fourier (ARSF) 
 * coefficients computed with downward recursion of modified Bessel functions.
 * Values are computed using LUT!
 */
void updateARSFCoeffRecursDownLUT(double lambda,double phi,int n,Eigen::VectorXd& coeffs);

/** Evaluates the Angular Radon Spectrum Fourier (ARSF) approximation according to "naive"
 * evaluation.
 */
double evaluateARSFRaw(const Eigen::VectorXd& coeffs,double theta);

/** Computes coeffients of ARSF Correlation. 
 */
void computeARSFCorr(const Eigen::VectorXd& arsfSrc,const Eigen::VectorXd& arsfDst,Eigen::VectorXd& arsfCor);

// --------------------------------------------------------
// OPTIMIZATION AND INTERVAL FUNCTIONS
// --------------------------------------------------------

/** Computes lower and upper bounds of cosine function on a given interval.
 */
void findCosLU(double a,double b,double& cmin,double& cmax);

/** Computes lower and upper bounds of an ARSF function (represented by its coefficients)
 * on a given interval.
 */
void findARSFLU(const Eigen::VectorXd& coeffs,double theta0,double theta1,double& arsfMin,double& arsfMax);

/** Computes the maximum of an ARSF function using Branch & Bound (BB) approach. 
 */
void findARSFMaxBB(const Eigen::VectorXd& coeffs,double theta0,double theta1,double thetaToll,double arsfToll,double& thetaMax,double& arsfMax);


// --------------------------------------------------------
// ANGULAR RADON SPECTRUM CLASS
// --------------------------------------------------------

/** Class for computing Angular Radon Spectum (ARS) of a set of points.
 */
class AngularRadonSpectrum
{
public:
  /** Default constructor. 
   */
  AngularRadonSpectrum() 
   : pairData_(), coeffsDirect_(), coeffsRecursDown_(), sigma_(1.0), pointNum_(0), arsfOrder_(0),
     thetaToll_(M_PI/180.0*0.5), threadNumOMP_(4), pnebiLut_()
  { }
 
  /** Sets the order of Fourier approximation of ARS.
   * WARNING: It resets the coefficients!
   */
  void setARSFOrder(int n)
  {
    arsfOrder_ = n;
    coeffsDirect_.resize(2*n+2);
    coeffsRecursDown_.resize(2*n+2);
  }

  /** Sets the standard deviation of isotropic gaussians used in point density
   * representation.
   */
  void setSigma(double sigma) { sigma_ = sigma; } 

  /** Sets the maximum tollerance on theta during the computation of maximum.
   */
  void setThetaToll(double thetaToll) { thetaToll_ = thetaToll; }

  /** Sets the number of thread used by OpenMP routines. 
   */
  void setThreadNumOMP(int tno) { threadNumOMP_ = tno; }

  /** Returns const reference to ARSF coefficients obtained from direct computation.
   */
  const Eigen::VectorXd& coeffsDirect() const { return coeffsDirect_; }

  /** Returns const reference to ARSF coefficients obtained from downward recursion. 
   */
  const Eigen::VectorXd& coeffsRecursDown() const { return coeffsRecursDown_; }

  /** Inserts the given points and computes all the data about point pairs.
   */
  template <typename PointIt>
  void insertPoints(PointIt beg,PointIt end);

  /** Initializes LUT (the LUT is used by initARSFRecursDownLUT).
   */
  void initLUT() 
  { 
    pnebiLut_.init(arsfOrder_,0.005);
  }

  /** Initializes the coefficients of Fourier series expansion of Radon Spectrum
   * arrested to order n-th. 
   * It internally computes coefficient vector (of size 2*n+2):
   *  coeffs_ = [ a_0  b_0 a_1  b_1 ... a_n b_n ]
   */
  void initARSFDirect();

  /** Initializes the coefficients of Fourier series expansion of Radon Spectrum
   * arrested to order n-th. 
   * It internally computes coefficient vector (of size 2*n+2):
   *  coeffs_ = [ a_0  b_0 a_1  b_1 ... a_n b_n ]
   */
  void initARSFRecursDown();

  /** Initializes the coefficients of Fourier series expansion of Radon Spectrum
   * arrested to order n-th. 
   * It internally computes coefficient vector (of size 2*n+2):
   *  coeffs_ = [ a_0  b_0 a_1  b_1 ... a_n b_n ]
   * It uses LUT to speed up computation.
   */
  void initARSFRecursDownLUT();

  /** Initializes and computes the ARS coefficients using a parallel routines  
   * of OpenMP (OMP)
   */
  void initARSFRecursDownOMP(const std::vector<Eigen::Vector2d>& points);

  /** Initializes the coefficients of Fourier series expansion of Radon Spectrum
   * arrested to order n-th. 
   * It does it without storing PairData. This function has been created to remove 
   * potential inefficiency. 
   */
  void initFromPoints(const std::vector<Eigen::Vector2d>& points);

  /** Direct evaluation of GMM-RS in exponential of cosines form:
   *   f(theta) = exp(-lambda - lambda * cos(theta - phi2))
   */
  double evalExpCos(double theta) const;

  /** Evaluates the ARSF using the coefficients obtained from direct computation.
   */
  double evalARSFDirect(double theta) const;

  /** Evaluates the ARSF using the coefficients obtained from downward recursion. 
   */
  double evalARSFRecursDown(double theta) const;

  /** Finds the maximum of ARSF.
   */
  double findMax() const 
  {
    double arsfMax, thetaMax;
    findARSFMaxBB(coeffsRecursDown_,0,M_PI,thetaToll_,10.0,thetaMax,arsfMax);
    return arsfMax;
  }

  /** Finds the maximum of ARSF.
   */
  double findMax(double& thetaMax) const 
  {
    double arsfMax;
    findARSFMaxBB(coeffsRecursDown_,0,M_PI,thetaToll_,10.0,thetaMax,arsfMax);
    return arsfMax;
  }

private:
  // Computes
  struct PairData
  {
    int i;
    int j;
    double lambda;
    double phi;
  };

  // Members
  std::vector<PairData> pairData_;
  Eigen::VectorXd coeffsDirect_;
  Eigen::VectorXd coeffsRecursDown_;
  double sigma_; 
  int pointNum_;
  int arsfOrder_;
  double thetaToll_;
  int threadNumOMP_;
  // PNEBI LUT
  PNEBILUT pnebiLut_;
};

template <typename PointIt>
void AngularRadonSpectrum::insertPoints(PointIt beg,PointIt end)
{
  PairData pdata;
  PointIt itNxt, endPrv;
  double x, y, factor;
  int i, j, n;

  if (beg == end) {
    std::cerr << __FILE__ << "," << __LINE__ << ": no points!" << std::endl;
    return;
  }
  factor = 1.0 / (8.0 * sigma_ * sigma_);

  // end1: position before the last one (if there is any!)
  endPrv = end;
  std::advance(endPrv,-1);
  // Visits each point
  i = 0;
  n = std::distance(beg,end);
  pairData_.clear();
  pairData_.reserve(n*(n-1));
  for (PointIt it1 = beg; it1 != endPrv; ++it1, ++i) {
    itNxt = it1;
    std::advance(itNxt,1);
    j = i+1;
    for (PointIt it2 = itNxt; it2 != end; ++it2, ++j) {
      pdata.i = i;
      pdata.j = j;
      x = it2->x() - it1->x();
      y = it2->y() - it1->y();
      pdata.lambda = factor * (x * x + y * y);
      pdata.phi = atan2(y,x);
      pairData_.push_back(pdata);
//      std::cout << "pair " << i << "," << j << ": lambda " << pdata.lambda << ", phi " << pdata.phi << std::endl;
    }
  }
  // Number of points
  pointNum_ = i + 1;
}

// ----------------------------------------------
// ARS CORRELATION CLASS
// ----------------------------------------------

/** Class for computing the best angle correlation of two Angular Radon Spectra (ARS).
 */
class ARSCorrelation
{
public:
  /** Constructor.
   */
  ARSCorrelation();

  /** Initializes the correlation of two ARS. 
   * The correlation function is defined as the integral:
   *
   *  corr(\delta) = \frac{1}{\pi} \int_{t=0}^{\pi} arsSrc(t+\delta) * arsDst(t) * dt
   */
  void init(const AngularRadonSpectrum& arsSrc,const AngularRadonSpectrum& arsDst);

  /** Another initialization of correlation using coefficients. 
   */
  void init(const Eigen::VectorXd& coeffsSrc,const Eigen::VectorXd& coeffsDst);

  /** Sets the tollerance on the estimation of angle with max correlation. 
   */
  void setAngleToll(double toll) { angleToll_ = toll; }

  /** Returns a reference to the coefficients of Fourier expanson on correlation. 
   */
  const Eigen::VectorXd& coeffs() const { return coeffsCorr_; }

  /** Gets the angle of maximum correlation (once computed!).
   */
  double getAngleMax() const { return angleMax_; }

  /** Gets the angle of maximum correlation (once computed!).
   */
  double getCorrMax() const { return corrMax_; }

  /** Estimates the angle that maximize correlation.
   */
  double findAngleMax();

  /** Evaluates the value of correlation based on Fourier approximation. 
   */
  double eval(double angle) const;

private:
  Eigen::VectorXd coeffsCorr_;
  double corrMax_;
  double angleMax_;
  double angleToll_;
};

}  // end of namespace

