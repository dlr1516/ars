/**
 * ARS - Angular Radon Spectrum 
 * Copyright (C) 2017-2020 Dario Lodi Rizzini.
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
#include <ars/MortonOrderedPoints.h>

namespace ars {
    
    template class MortonOrderedPoints<2, 4, float>;
    template class MortonOrderedPoints<2, 4, double>;
    template class MortonOrderedPoints<2, 8, float>;
    template class MortonOrderedPoints<2, 8, double>;
    template class MortonOrderedPoints<2, 16, float>;
    template class MortonOrderedPoints<2, 16, double>;
    template class MortonOrderedPoints<2, 24, float>;
    template class MortonOrderedPoints<2, 24, double>;
    template class MortonOrderedPoints<2, 32, float>;
    template class MortonOrderedPoints<2, 32, double>;
    
    template class MortonOrderedPoints<3, 4, float>;
    template class MortonOrderedPoints<3, 4, double>;
    template class MortonOrderedPoints<3, 8, float>;
    template class MortonOrderedPoints<3, 8, double>;
    template class MortonOrderedPoints<3, 16, float>;
    template class MortonOrderedPoints<3, 16, double>;
    template class MortonOrderedPoints<3, 24, float>;
    template class MortonOrderedPoints<3, 24, double>;
    template class MortonOrderedPoints<3, 32, float>;
    template class MortonOrderedPoints<3, 32, double>;
};
