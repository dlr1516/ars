#ifndef ARS_KERNEL_H
#define ARS_KERNEL_H

#include <iostream>
#include <ars/definitions.h>

namespace ars {

/** 
 * Class functor for the evaluation of an ARS kernel function derived from a pair 
 * of Gaussian distributions. 
 * The ARS kernel function has the form:
 * 
 * 1 / sqrt() 
 */
class ArsKernel {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /**
   * Constructor of an ARS kernel 
   */
  ArsKernel(const Vector2& mean1,const Matrix2& covar1,const Vector2& mean2,const Matrix2& covar2);

  /**
   * Constructor of an ARS kernel 
   */
  ArsKernel(double lamda, double sigma, double phaseNum, double phaseDen);

  /**
   */
  double evaluate(double theta) const;

private:
  double delta_;
  double sigma_;
  double phaseNum_;
  double phaseDen_;
};

} // end of namespace

#endif 

