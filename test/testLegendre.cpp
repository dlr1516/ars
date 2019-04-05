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
#include <iostream>
#include <ars/functions.h>
#include <ars/ars2d.h>
#include <ars/thirdparty/gnuplot-iostream.h>
#include <ars/Profiler.h>
#include <ars/ParamMap.h>

int main(int argc, char** argv) {
    ars::ParamMap paramMap;
    
    
    paramMap.read(argc, argv);
    
    return 0;
}

