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

namespace ars {

typedef Eigen::Vector2d Point2;

typedef std::vector<Point2, Eigen::aligned_allocator<Point2> > Point2Vector;

} // end of namespace

#endif /* DEFINITIONS_H */

