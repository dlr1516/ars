#include <emotion/AngularRadonSpectrum.h>
#include <boost/numeric/interval.hpp>
#include <queue>

#include <omp.h>

#define STREAM_BOUND_INTERVAL(X)  "[" << (180.0/M_PI*(X).theta0) << "," << (180.0/M_PI*(X).theta1) << "] x [" << (X).lower << "," << (X).upper << "]"

namespace emotion
{

// --------------------------------------------------------
// ARSF FUNCTIONS
// --------------------------------------------------------

void updateARSFCoeffDirect(double lambda,double phi,int n,Eigen::VectorXd& coeffs)
{
  double pnebi, sgn;

  // Fourier Coefficients 
  if (coeffs.size() != 2*n + 2) {
    std::cerr << __FILE__ << "," << __LINE__ << ": invalid size of Fourier coefficients vector " << coeffs.size() << " should be " << (2*n+2) << std::endl;
    return;
  }

  // Computes coefficients: the new values are added or set according to flag
  coeffs(0) += 0.5 * evaluatePNEBI(0,lambda);
  sgn = -1.0;
  for (int  k = 1; k <= n; ++k) {
    pnebi = evaluatePNEBI(k,lambda);
    //std::cout << "k " << k << ", pnebi " << pnebi << " coeffs(" << (2*k) << ") and coeffs(" << (2*k+1) << "), size " << coeffs.size() << std::endl;
    coeffs(2*k)   += pnebi * sgn * cos(2.0 * k * phi);
    coeffs(2*k+1) += pnebi * sgn * sin(2.0 * k * phi);
    sgn = -sgn;
  }
}

void updateARSFCoeffRecursUp(double lambda,double phi,int n,Eigen::VectorXd& coeffs)
{
  double cphi2, sphi2, cphi4, sphi4;
  double ak, bk, ak1, bk1, ak2, bk2;
  double modk, modk1, multK1;

  // Fourier Coefficients 
  if (coeffs.size() != 2*n + 2) {
    std::cerr << __FILE__ << "," << __LINE__ << ": invalid size of Fourier coefficients vector " << coeffs.size() << " should be " << (2*n+2) << std::endl;
    return;
  }

  cphi2 = cos(2.0 * phi);
  sphi2 = sin(2.0 * phi);
  cphi4 = cos(4.0 * phi);
  sphi4 = sin(4.0 * phi);
  ak2 = evaluatePNEBI(0,lambda);
  bk2 = 0.0;
  ak1 = -evaluatePNEBI(1,lambda) * cphi2;
  bk1 = -evaluatePNEBI(1,lambda) * sphi2;
  coeffs(0) += 0.5 * ak2;
  coeffs(1) += bk2;
  coeffs(2) += ak1;
  coeffs(3) += bk1;
  modk1 = ak1*ak1 + bk1*bk1;
  for (int k = 2; k <= n; ++k) {
    multK1 = 2.0 * (k-1) / lambda;
    ak = ak2 * cphi4 - bk2 * sphi4 + multK1 * (ak1 * cphi2 - bk1 * sphi2);
    bk = ak2 * sphi4 + bk2 * cphi4 + multK1 * (ak1 * sphi2 + bk1 * cphi2);
    modk = ak*ak + bk*bk;
    coeffs(2*k  ) += ak;
    coeffs(2*k+1) += bk;
    std::swap(ak2,ak1);
    std::swap(bk2,bk1);
    std::swap(ak1,ak);
    std::swap(bk1,bk);
  }
}

/** Computes the contribution of term (lambda,phi) to the Fourier coefficients of ARS
 * computed with downward recursion of modified Bessel functions.
 */
void updateARSFCoeffRecursDown(double lambda,double phi,int n,Eigen::VectorXd& coeffs)
{
  Eigen::VectorXd pnebis(n+1);
  double sgn, cth2, sth2, cth, sth, ctmp, stmp;

  // Fourier Coefficients 
  if (coeffs.size() != 2*n + 2) {
    std::cerr << __FILE__ << "," << __LINE__ << ": invalid size of Fourier coefficients vector " << coeffs.size() << " should be " << (2*n+2) << std::endl;
    return;
  }

  evaluatePNEBIVector(n,lambda,pnebis);
  coeffs(0) += 0.5 * pnebis(0);
  sgn = -1.0;
  cth2 = cos(2.0*phi);
  sth2 = sin(2.0*phi);
  cth = cth2;
  sth = sth2;
  for (int  k = 1; k <= n; ++k) {
//    coeffs(2*k)   += pnebis(k) * sgn * cos(2.0 * k * phi);
//    coeffs(2*k+1) += pnebis(k) * sgn * sin(2.0 * k * phi);
    coeffs(2*k)   += pnebis(k) * sgn * cth;
    coeffs(2*k+1) += pnebis(k) * sgn * sth;
    sgn = -sgn;
    ctmp = cth2 * cth - sth2 * sth;
    stmp = sth2 * cth + cth2 * sth;
    cth = ctmp;
    sth = stmp;
  }
}

void updateARSFCoeffRecursDownLUT(double lambda,double phi,int n,const PNEBILUT& pnebiLUT,Eigen::VectorXd& coeffs)
{
  Eigen::VectorXd pnebis(n+1);
  double sgn, cth2, sth2, cth, sth, ctmp, stmp;

  // Fourier Coefficients 
  if (coeffs.size() != 2*n + 2 || pnebiLUT.getOrderMax() < n) {
    std::cerr << __FILE__ << "," << __LINE__ << ": one of these conditions failed:"
     << "\n  size of Fourier coefficients vector " << coeffs.size() << " should be " << (2*n+2)
     << "\n  LUT max order is " << pnebiLUT.getOrderMax() << " >= " << n 
     << std::endl;
    return;
  }
  
  pnebiLUT.eval(lambda,pnebis);
  coeffs(0) += 0.5 * pnebis(0);
  sgn = -1.0;
  //cth2 = cos(2.0*phi);
  //sth2 = sin(2.0*phi);
  // Fast approximation of cos() and sin()
  emotion::fastCosSin(2.0*phi,cth2,sth2);
  cth = cth2;
  sth = sth2;
  for (int  k = 1; k <= n; ++k) {
    coeffs(2*k)   += pnebis(k) * sgn * cth;
    coeffs(2*k+1) += pnebis(k) * sgn * sth;
    sgn = -sgn;
    ctmp = cth2 * cth - sth2 * sth;
    stmp = sth2 * cth + cth2 * sth;
    cth = ctmp;
    sth = stmp;
  }
}


double evaluateARSFRaw(const Eigen::VectorXd& coeffs,double theta)
{
  double val, cth2, sth2, cth, sth, ctmp, stmp;
  int n;

  if (coeffs.size() % 2 != 0) {
    std::cerr << __FILE__ << "," << __LINE__ << ": the number of coefficients must be even: found " << coeffs.size() << std::endl;
  }
  n = (coeffs.size() / 2) - 1;

  cth2 = cos(2.0*theta);
  sth2 = sin(2.0*theta);
  cth = 1.0;
  sth = 0.0;
  val = 0.0;

  for (int k = 0; k <= n; ++k) {
    val += coeffs(2*k) * cth + coeffs(2*k+1) * sth;
    ctmp = cth2 * cth - sth2 * sth;
    stmp = sth2 * cth + cth2 * sth;
    cth = ctmp;
    sth = stmp;
  }
  return val;
}

void computeARSFCorr(const Eigen::VectorXd& arsfSrc,const Eigen::VectorXd& arsfDst,Eigen::VectorXd& arsfCor)
{
  int n;

  if (arsfSrc.size() % 2 != 0 || arsfDst.size() % 2 != 0 || arsfSrc.size() != arsfDst.size()) {
    std::cerr << __FILE__ << "," << __LINE__ << ": the number of ARSF coefficients must be even and equal: arsfSrc " << arsfSrc.size() 
      << ", arsfDst " << arsfDst.size() << std::endl;
    return;
  }
  n = (arsfSrc.size() / 2) - 1;  

  arsfCor.resize(2*n+2);
  for (int k = 0; k <= n; ++k) {
    arsfCor(2*k)   = 0.5 * (arsfSrc(2*k) * arsfDst(2*k)   + arsfSrc(2*k+1) * arsfDst(2*k+1));
    arsfCor(2*k+1) = 0.5 * (arsfSrc(2*k) * arsfDst(2*k+1) - arsfSrc(2*k+1) * arsfDst(2*k));
  }
}

// --------------------------------------------------------
// OPTIMIZATION AND INTERVAL FUNCTIONS
// --------------------------------------------------------

void findCosLU(double a,double b,double& cmin,double& cmax)
{
  double amod, bmod;

  if (a > b) std::swap(a,b);

  if (b - a >= 2.0 * M_PI) {
    cmin = -1.0;
    cmax = +1.0;
  }
  else {
    // Normalizes to circular interval [0, 2*M_PI[
    amod = fmod(a,2.0*M_PI);
    if (amod < 0.0) amod += 2.0*M_PI;
    bmod = fmod(b,2.0*M_PI);
    if (bmod < 0.0) bmod += 2.0*M_PI;
    // Case bmod < amod: for example [300,30[ deg: angle 0 is included.
    if (bmod < amod) {
      cmax = +1.0;
      if (bmod < M_PI && M_PI < amod) {
        cmin = std::min(cos(amod),cos(bmod));
      }
      else {
        cmin = -1.0;
      }
//      if (M_PI < bmod || amod < M_PI) {
//        cmin = -1.0;
//      }
//      else {
//        cmin = std::min(cos(amod),cos(bmod));
//      }
    }
    else {
      cmax = std::max(cos(amod),cos(bmod));
      if (amod < M_PI && M_PI < bmod) {
        cmin = -1.0;
      }
      else {
        cmin = std::min(cos(amod),cos(bmod));
      }
    }
  }
}

void findARSFLU(const Eigen::VectorXd& coeffs,double theta0,double theta1,double& arsfMin,double& arsfMax)
{
  double amplitude, phase, sinusoidMin, sinusoidMax;
  int n, i0, i1;

  if (coeffs.size() % 2 != 0) {
    std::cerr << __FILE__ << "," << __LINE__ << ": the number of coefficients must be even: found " << coeffs.size() << std::endl;
  }
  n = (coeffs.size() / 2) - 1;

  if (theta1 < theta0) {
    std::cerr << __FILE__ << "," << __LINE__ << ": invalid interval [" << theta0 << "," << theta1 << "]: swapping endpoints to continue" << std::endl;
    std::swap(theta0,theta1);
  }

  // arsfMin and arsfMax initialized with constant component
  arsfMin = coeffs(0);
  arsfMax = coeffs(0);
  for (int k = 1; k <= n; ++k) {
    // t_k = a_k * cos(2*k*theta) + b_k * sin(2*k*theta) = amplitude * cos(2*k*theta - phase)
    // Period of k-th terms is M_PI / k. 
    amplitude = sqrt(coeffs(2*k) * coeffs(2*k) + coeffs(2*k+1) * coeffs(2*k+1));
    phase = atan2(coeffs(2*k+1),coeffs(2*k));
    //std::cout << "k " << k << ", amplitude " << amplitude << ", phase[deg] " << (180.0/M_PI*phase) << std::endl;
    // If the [theta0,theta1] interval is larger than period, then the whole sinusoid amplitude is considered.
    // Otherwise, a more refined evaluation is performed.
    findCosLU(2*k*theta0 - phase,2*k*theta1 - phase,sinusoidMin,sinusoidMax);
    arsfMin += amplitude * sinusoidMin;
    arsfMax += amplitude * sinusoidMax;
  }
}

void findARSFMaxBB(const Eigen::VectorXd& coeffs,double theta0,double theta1,double thetaToll,double arsfToll,double& thetaMax,double& arsfMax)
{
  struct IntervalBound
  {
    double theta0;
    double theta1;
    double lower;
    double upper;
  };
  struct IntervalBoundCmp
  {
    bool operator()(const IntervalBound& ib0,const IntervalBound& ib1) const { return (ib0.upper < ib1.upper); }
  } cmp;
  std::priority_queue<IntervalBound,std::vector<IntervalBound>,IntervalBoundCmp> queue;
  IntervalBound curr, left, right, global;
  
  global.theta0 = theta0;
  global.theta1 = theta1;
  findARSFLU(coeffs,global.theta0,global.theta1,global.lower,global.upper);
  queue.push(global);
  //std::cout << "init " << STREAM_BOUND_INTERVAL(global) << std::endl;

  while (!queue.empty()) {
    curr = queue.top();
    queue.pop();
    // It processes the interval curr only if the further analysis can improve the current solution. 
    // In practice, if curr.upper < global.lower, then curr is ignored.
    if (curr.upper >= global.lower) {
      // Updates: the current candidate interval to contain global maximum with the lower and upper bounds
      // of global maximum
      if (global.lower <= curr.lower) {
        global.theta0 = curr.theta0;
        global.theta1 = curr.theta1;
        global.lower = curr.lower;
        global.upper = curr.upper;
      }
//      std::cout << "curr " << STREAM_BOUND_INTERVAL(curr) << std::endl;
//      std::cout << "  global " << STREAM_BOUND_INTERVAL(global) << std::endl; 
      // Splits curr into intervals left and right and computes bounds on both of them
      if (curr.theta1 - curr.theta0 > thetaToll) {
        left.theta0 = curr.theta0;
        left.theta1 = 0.5 * (curr.theta0 + curr.theta1);
        findARSFLU(coeffs,left.theta0,left.theta1,left.lower,left.upper);
        right.theta0 = left.theta1;
        right.theta1 = curr.theta1;
        findARSFLU(coeffs,right.theta0,right.theta1,right.lower,right.upper);  
        queue.push(left);
        queue.push(right);
      }
    }
    else {
//      std::cout << "  discarding " << STREAM_BOUND_INTERVAL(curr) << std::endl;
    }
  }
  thetaMax = 0.5 * (global.theta0 + global.theta1);
  arsfMax = 0.5 * (global.lower + global.upper);
}

// ----------------------------------------------
// ANGULAR RADON SPECTRUM CLASS
// ----------------------------------------------

void AngularRadonSpectrum::initARSFDirect()
{
  coeffsDirect_.fill(0.0);
  //std::cout << __FILE__ << "," << __LINE__ << ": \n" << coeffsDirect_.transpose() << std::endl;
//  std::cout << __FILE__ << "," << __LINE__ << ": coeffsDirect_\n" << std::endl;
  for (auto& p : pairData_) {
    updateARSFCoeffDirect(p.lambda,p.phi,arsfOrder_,coeffsDirect_);
//    std::cout << coeffsDirect_.transpose() << std::endl;
  }
  coeffsDirect_(0) += 0.5 * pointNum_;
  coeffsDirect_ = coeffsDirect_ / (sigma_ * sqrt(M_PI));
}

void AngularRadonSpectrum::initARSFRecursDown()
{
  coeffsRecursDown_.fill(0.0);
  //std::cout << __FILE__ << "," << __LINE__ << ": \n" << coeffsRecursDown_.transpose() << std::endl;
//  std::cout << __FILE__ << "," << __LINE__ << ": coeffsRecursDown_\n" << std::endl;
  for (auto& p : pairData_) {
    updateARSFCoeffRecursDown(p.lambda,p.phi,arsfOrder_,coeffsRecursDown_);
//    std::cout << coeffsRecursDown_.transpose() << std::endl;
  }
  coeffsRecursDown_(0) += 0.5 * pointNum_;
  coeffsRecursDown_ = coeffsRecursDown_ / (sigma_ * sqrt(M_PI));
}

void AngularRadonSpectrum::initARSFRecursDownLUT()
{
  // Check PNEBI function LUT
  if (pnebiLut_.getOrderMax() < arsfOrder_) {
    std::cerr << __FILE__ << "," << __LINE__ << ": LUT not initialized to right order. Initialized now." << std::endl;
    pnebiLut_.init(arsfOrder_,0.005);
  }

  // Checking Fourier coefficients vector
  if (coeffsRecursDown_.size() != 2*arsfOrder_ + 2) {
    std::cerr << __FILE__ << "," << __LINE__ << ": size of Fourier coefficients vector " << coeffsRecursDown_.size() 
      << " should be " << (2*arsfOrder_+2) << ": resized" << std::endl;
    coeffsRecursDown_.resize(2*arsfOrder_+2);
  }

  // Variables
  Eigen::VectorXd pnebis(arsfOrder_+1);
  double sgn, cth2, sth2, cth, sth, ctmp, stmp;

  coeffsRecursDown_.fill(0.0);
  for (auto& p : pairData_) {
    // Copies values of pnebiLUT into pnebis vector
    pnebiLut_.eval(p.lambda,pnebis);
    coeffsRecursDown_(0) += 0.5 * pnebis(0);
    sgn = -1.0;
    emotion::fastCosSin(2.0*p.phi,cth2,sth2);
    cth = cth2;
    sth = sth2;
    // Computing coefficients of all orders
    for (int  k = 1; k <= arsfOrder_; ++k) {
      coeffsRecursDown_(2*k)   += pnebis(k) * sgn * cth;
      coeffsRecursDown_(2*k+1) += pnebis(k) * sgn * sth;
      sgn = -sgn;
      ctmp = cth2 * cth - sth2 * sth;
      stmp = sth2 * cth + cth2 * sth;
      cth = ctmp;
      sth = stmp;
    }
    // Hard-coded update of coefficients avoids unnecessary memory allocatation/deallocation
    // inside the following function:
    // updateARSFCoeffRecursDownLUT(p.lambda,p.phi,arsfOrder_,pnebiLut_,coeffsRecursDown_);
  }
  coeffsRecursDown_(0) += 0.5 * pointNum_;
  coeffsRecursDown_ = coeffsRecursDown_ / (sigma_ * sqrt(M_PI));
  
//  coeffsRecursDown_.fill(0.0);
//  for (auto& p : pairData_) {
//    updateARSFCoeffRecursDownLUT(p.lambda,p.phi,arsfOrder_,pnebiLut_,coeffsRecursDown_);
//  }
//  coeffsRecursDown_(0) += 0.5 * pointNum_;
//  coeffsRecursDown_ = coeffsRecursDown_ / (sigma_ * sqrt(M_PI));
}

void AngularRadonSpectrum::initARSFRecursDownOMP(const std::vector<Eigen::Vector2d>& points)
{
  double factor;
  int blockSize;
  // OMP does not allows class member as shared variables: 
  // http://stackoverflow.com/questions/4610656/why-is-class-member-variable-x-not-allowed-to-be-sharedx-in-openmp
  // A copy of the variables is required!
  Eigen::VectorXd coeffsGlobal(coeffsRecursDown_.size());
  int pointNumCopy = points.size();

  pointNum_ = points.size();
  blockSize = pointNum_ / threadNumOMP_ + ((pointNum_ % threadNumOMP_ != 0));
  factor = 1.0 / (8.0 * sigma_ * sigma_);
  coeffsGlobal.fill(0.0);

  // http://bisqwit.iki.fi/story/howto/openmp/
  //#pragma omp parallel num_threads(threadNumOMP_) shared(points,coeffsRecursDown_,pointNum_,arsfOrder_) private(factor,blockSize) 
//  {
//    //#pragma omp for ordered schedule(dynamic) 
//    #pragma omp parallel for ordered schedule(dynamic) num_threads(threadNumOMP_) shared(points,factor,blockSize,coeffsGlobal,pointNumCopy) 
//    for (int b = 0; b < threadNumOMP_; ++b) {
//      Eigen::VectorXd coeffsLocal(coeffsGlobal.size());
//      coeffsLocal.fill(0.0);
//      double x, y, lambda, phi;
//      int idxBeg = b * blockSize;
//      int idxEnd = std::min((b+1)*blockSize, pointNumCopy-1);
//      for (int i = idxBeg; i < idxEnd; ++i) {
//        for (int j = i+1; j < pointNumCopy; ++j) {
//          x = points[i].x() - points[j].x();
//          y = points[i].y() - points[j].y();
//          lambda = factor * (x * x + y * y);
//          phi = atan2(y,x);
//          updateARSFCoeffRecursDown(lambda,phi,arsfOrder_,coeffsLocal);
//        }
//      }
//      for (int i = 0; i < coeffsGlobal.size(); ++i) {
//        #pragma omp atomic
//        coeffsGlobal(i) += coeffsLocal(i);
//      }
//    }
//  }

  #pragma omp parallel num_threads(threadNumOMP_) shared(points,factor,blockSize,coeffsGlobal,pointNumCopy) 
  {
    Eigen::VectorXd coeffsLocal(coeffsGlobal.size());
    coeffsLocal.fill(0.0);
    double x, y, lambda, phi;
    #pragma omp for ordered schedule(static) 
    for (int i = 0; i < pointNumCopy; ++i) {
      for (int j = i+1; j < pointNumCopy; ++j) {
        x = points[i].x() - points[j].x();
        y = points[i].y() - points[j].y();
        lambda = factor * (x * x + y * y);
        phi = atan2(y,x);
        updateARSFCoeffRecursDown(lambda,phi,arsfOrder_,coeffsLocal);
      }
    }
    for (int i = 0; i < coeffsGlobal.size(); ++i) {
      #pragma omp atomic
      coeffsGlobal(i) += coeffsLocal(i);
    }
  }


  coeffsRecursDown_ = coeffsGlobal;
  coeffsRecursDown_(0) += 0.5 * pointNum_;
  coeffsRecursDown_ = coeffsRecursDown_ / (sigma_ * sqrt(M_PI));
}

void AngularRadonSpectrum::initFromPoints(const std::vector<Eigen::Vector2d>& points)
{
  double x, y, phi, lambda, factor;
  int n = points.size();

  factor = 1.0 / (8.0 * sigma_ * sigma_);
  for (int i = 0; i < n-1; ++i) {
    for (int j = i+1; j < n; ++j) {
      
    }
  }
}

// Direct evaluation of exponential of cosines
double AngularRadonSpectrum::evalExpCos(double theta) const
{
  double val = 0.0;
  double arg;
  for (auto& p : pairData_) {
    val += exp(-p.lambda * (1.0 + cos(2.0*(theta - p.phi))) );
  }
  val = (val + 0.5 * pointNum_) / (sigma_ * sqrt(M_PI));
  return val;
}

double AngularRadonSpectrum::evalARSFDirect(double theta) const
{
  return evaluateARSFRaw(coeffsDirect_,theta);
}

double AngularRadonSpectrum::evalARSFRecursDown(double theta) const
{
  return evaluateARSFRaw(coeffsRecursDown_,theta);
}

// ----------------------------------------------
// ARS CORRELATION CLASS
// ----------------------------------------------

ARSCorrelation::ARSCorrelation() : coeffsCorr_(), angleMax_(0.0), angleToll_(M_PI/180.0*1.0)
{
}

void ARSCorrelation::init(const Eigen::VectorXd& coeffsSrc,const Eigen::VectorXd& coeffsDst)
{
  computeARSFCorr(coeffsSrc,coeffsDst,coeffsCorr_);
}

void ARSCorrelation::init(const AngularRadonSpectrum& arsSrc,const AngularRadonSpectrum& arsDst)
{
  computeARSFCorr(arsSrc.coeffsRecursDown(),arsDst.coeffsRecursDown(),coeffsCorr_);
}

double ARSCorrelation::findAngleMax()
{
  findARSFMaxBB(coeffsCorr_,0,M_PI,angleToll_,10.0,angleMax_,corrMax_);
  return angleMax_;
}

double ARSCorrelation::eval(double angle) const
{
  return evaluateARSFRaw(coeffsCorr_,angle);
}

}  // end of namespace

