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
#ifndef ARS_DEFINITIONS_H
#define ARS_DEFINITIONS_H

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/StdVector>

#define ARS_PRINT(MSG) std::cout << __FILE__ << "," << __LINE__ << ": " << MSG << std::endl;

#define ARS_ERROR(MSG) std::cerr << __FILE__ << "," << __LINE__ << ": " << MSG << std::endl;

#define ARS_VAR1(X1) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) << std::endl;

#define ARS_VAR2(X1,X2) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) \
   << ", " << (#X2) << " " << (X2) << std::endl;

#define ARS_VAR3(X1,X2,X3) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) \
   << ", " << (#X2) << " " << (X2) << ", " << (#X3) << " " << (X3) << std::endl;

#define ARS_VAR4(X1,X2,X3,X4) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) \
   << ", " << (#X2) << " " << (X2) << ", " << (#X3) << " " << (X3) << ", " << (#X4) << " " << (X4) << std::endl;

#define ARS_VAR5(X1,X2,X3,X4,X5) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) \
   << ", " << (#X2) << " " << (X2) << ", " << (#X3) << " " << (X3) << ", " << (#X4) << " " << (X4) \
   << ", " << (#X5) << " " << (X5) << std::endl;

#define ARS_VAR6(X1,X2,X3,X4,X5,X6) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) \
   << ", " << (#X2) << " " << (X2) << ", " << (#X3) << " " << (X3) << ", " << (#X4) << " " << (X4) \
   << ", " << (#X5) << " " << (X5) << ", " << (#X6) << " " << (X6) << std::endl;

#define ARS_VAR7(X1,X2,X3,X4,X5,X6,X7) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X1) << " " << (X1) \
   << ", " << (#X2) << " " << (X2) << ", " << (#X3) << " " << (X3) << ", " << (#X4) << " " << (X4) \
   << ", " << (#X5) << " " << (X5) << ", " << (#X6) << " " << (X6) << ", " << (#X7) << " " << (X7) << std::endl;

#define ARS_ASSERT(COND) if (!(COND)) { std::cerr << __FILE__ << "," << __LINE__ << ": assertion failed on " << #COND << std::endl; exit(-1); }

#define ARS_ASSERT_VAR1(COND,X1) \
   if (!(COND)) { std::cerr << __FILE__ << "," << __LINE__ << ": assertion failed on " << #COND << std::endl; \
   ARS_VAR1(X1); exit(-1); }

#define ARS_ASSERT_VAR2(COND,X1,X2) \
   if (!(COND)) { std::cerr << __FILE__ << "," << __LINE__ << ": assertion failed on " << #COND << std::endl; \
   ARS_VAR2(X1,X2); exit(-1); }

#define ARS_ASSERT_VAR3(COND,X1,X2,X3) \
   if (!(COND)) { std::cerr << __FILE__ << "," << __LINE__ << ": assertion failed on " << #COND << std::endl; \
   ARS_VAR3(X1,X2,X3); exit(-1); }

#define ARS_ASSERT_VAR4(COND,X1,X2,X3,X4) \
   if (!(COND)) { std::cerr << __FILE__ << "," << __LINE__ << ": assertion failed on " << #COND << std::endl; \
   ARS_VAR4(X1,X2,X3,X4); exit(-1); }

#define ARS_ASSERT_VAR5(COND,X1,X2,X3,X4,X5) \
   if (!(COND)) { std::cerr << __FILE__ << "," << __LINE__ << ": assertion failed on " << #COND << std::endl; \
   ARS_VAR5(X1,X2,X3,X4,X5); exit(-1); }

#define ARS_ASSERT_VAR6(COND,X1,X2,X3,X4,X5,X6) \
   if (!(COND)) { std::cerr << __FILE__ << "," << __LINE__ << ": assertion failed on " << #COND << std::endl; \
   ARS_VAR6(X1,X2,X3,X4,X5,X6); exit(-1); }

#define ARS_ASSERT_VAR7(COND,X1,X2,X3,X4,X5,X6,X7) \
   if (!(COND)) { std::cerr << __FILE__ << "," << __LINE__ << ": assertion failed on " << #COND << std::endl; \
   ARS_VAR7(X1,X2,X3,X4,X5,X6,X7); exit(-1); }

namespace ars {

#if __cplusplus < 201703L

    static_assert(__cplusplus >= 201703L, "Versions before C++-17 are deprecated: add compiling option \"-DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_FLAGS='-std=c++17'\"");
    
    using Vector2 = Eigen::Vector2d;

    using VectorVector2 = std::vector<Vector2, Eigen::aligned_allocator<Vector2> >;

    using Matrix2 = Eigen::Matrix2d;
    
    using VectorMatrix2 = std::vector<Matrix2, Eigen::aligned_allocator<Matrix2> >;
    
#else 
    
    using Vector2 = Eigen::Vector2d;

    using VectorVector2 = std::vector<Vector2>;

    using Matrix2 = Eigen::Matrix2d;
    
    using VectorMatrix2 = std::vector<Matrix2>;
    
#endif

} // end of namespace

#endif /* ARS_DEFINITIONS_H */

