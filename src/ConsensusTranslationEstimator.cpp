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
#include <ars/ConsensusTranslationEstimator.h>

namespace ars {
    
    template <> class ConsensusTranslationEstimator<2, float>;
    template <> class ConsensusTranslationEstimator<2, double>;
    template <> class ConsensusTranslationEstimator<3, float>;
    template <> class ConsensusTranslationEstimator<3, double>;
    
} // end of namespace
