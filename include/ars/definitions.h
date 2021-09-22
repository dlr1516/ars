/**
 * ARS - Angular Radon Spectrum 
 * Copyright (C) 2017 Dario Lodi Rizzini.
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
#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/StdVector>

#define ARS_PRINT(MSG) std::cout << __FILE__ << "," << __LINE__ << ": " << MSG << std::endl;

#define ARS_ERROR(MSG) std::cerr << __FILE__ << "," << __LINE__ << ": " << MSG << std::endl;

#define ARS_VARIABLE(X1) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) << std::endl;

#define ARS_VARIABLE2(X1,X2) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) \
   << ", " << (#X2) << " " << (X2) << std::endl;

#define ARS_VARIABLE3(X1,X2,X3) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) \
   << ", " << (#X2) << " " << (X2) << ", " << (#X3) << " " << (X3) << std::endl;

#define ARS_VARIABLE4(X1,X2,X3,X4) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) \
   << ", " << (#X2) << " " << (X2) << ", " << (#X3) << " " << (X3) << ", " << (#X4) << " " << (X4) << std::endl;

#define ARS_VARIABLE5(X1,X2,X3,X4,X5) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) \
   << ", " << (#X2) << " " << (X2) << ", " << (#X3) << " " << (X3) << ", " << (#X4) << " " << (X4) \
   << ", " << (#X5) << " " << (X5) << std::endl;

#define ARS_VARIABLE6(X1,X2,X3,X4,X5,X6) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) \
   << ", " << (#X2) << " " << (X2) << ", " << (#X3) << " " << (X3) << ", " << (#X4) << " " << (X4) \
   << ", " << (#X5) << " " << (X5) << ", " << (#X6) << " " << (X6) << std::endl;

#define ARS_VARIABLE7(X1,X2,X3,X4,X5,X6,X7) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) \
   << ", " << (#X2) << " " << (X2) << ", " << (#X3) << " " << (X3) << ", " << (#X4) << " " << (X4) \
   << ", " << (#X5) << " " << (X5) << ", " << (#X6) << " " << (X6) << ", " << (#X7) << " " << (X7) << std::endl;

#define ARS_ASSERT(COND) if (!(COND)) { std::cerr << __FILE__ << "," << __LINE__ << ": assertion failed on " << #COND << std::endl; exit(-1); }


namespace ars {

    typedef Eigen::Vector2d Vector2;

    typedef std::vector<Vector2, Eigen::aligned_allocator<Vector2> > VectorVector2;

    typedef Eigen::Matrix2d Matrix2;
    
    typedef std::vector<Matrix2, Eigen::aligned_allocator<Matrix2> > VectorMatrix2;

} // end of namespace

#endif /* DEFINITIONS_H */

