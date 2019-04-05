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

#include <Eigen/Dense>
#include <Eigen/StdVector>

#define ARS_PRINT(X) std::cout << __FILE__ << "," << __LINE__ << ": " << (#X) << " " << (X) << std::endl;

#define ARS_ASSERT(X) \
  if (!(X)) { \
    std::cerr << __FILE__ << "," << __LINE__ << ": assertion failed " << #X << " " << (X) << std::endl; \
    exit(-1); \
  } 

namespace ars {

    typedef Eigen::Vector2d Vector2;

    typedef std::vector<Vector2, Eigen::aligned_allocator<Vector2> > Vector2Vector;

    typedef Eigen::Matrix2d Matrix2;

    typedef Eigen::Vector3d Vector3;

    typedef std::vector<Vector3, Eigen::aligned_allocator<Vector3> > Vector3Vector;

    typedef Eigen::Matrix3d Matrix3;

} // end of namespace

#endif /* DEFINITIONS_H */

