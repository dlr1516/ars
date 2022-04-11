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
#ifndef ARS_FUNCTIONS_H
#define ARS_FUNCTIONS_H

#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
//#include <Eigen/Dense>

namespace ars {

    const double PNEBI_ARG_MAX = 600.0;
    const double BIG_NUM = 1.0e+10;
    const double SMALL_NUM = 1.0e-10;

    // --------------------------------------------------------
    // COS-SIN FAST EVALUATION
    // --------------------------------------------------------

    /** Computes sine and cosine values using a parabolic approximation of sine.
     *    y = 16.0 * xn * (abs(xn) - 0.5)
     * where xn = x/(2*M_PI) is the normalized value of x over the period. 
     *
     * Code suggested here:
     *  http://forum.devmaster.net/t/fast-and-accurate-sine-cosine/9648/6
     *  http://stackoverflow.com/questions/18662261/fastest-implementation-of-sine-cosine-and-square-root-in-c-doesnt-need-to-b
     */
    void fastCosSin(double x, double& c, double& s);
    
    /**
     * Computes atan() using a polynomial approximation on interval [-1,1]. See:
     * 
     *  Abramowitz, Stegun, "Handbook of Mathematical Functions", 1965
     * 
     * @param x the argument that must be in interval [-1.0, 1.0]
     * @return the value of atan
     */
    double fastAtan(double x);
    
    /**
     * Computes atan2() using fastAtan(). It uses clever interval as suggested by:
     *    https://www.dsprelated.com/showarticle/1052.php
     * (one comment in particular improved the discussed code). 
     * @param x
     * @param y
     * @return the approximate atan2. 
     */
    double fastAtan2(double x, double y);
    
    /**
     * Computes the value of the given Fourier series at the given point theta. 
     * The user must provide the vector of serie coefficients:
     *   S(x) = \sum_{i=0}^{n} ( coeffs[2*i] * cos(i*x) + coeffs[2*i+1] * sin(i*x) )
     * @param coeffs the vector of coefficiens (vector size must be an even number!)
     * @param theta the point where to compute the Fourier serie value
     * @return the value of the function
     */
    double evaluateFourier(const std::vector<double>& coeffs, double theta);

    // --------------------------------------------------------
    // PNEBI FUNCTIONS
    // PNEBI stands for Product of Negative Exponential and Bessel I, which is defined as
    // 
    //    PNEBI(k,x) = 2.0 * exp(-x) * besseli(k,x)
    // 
    // where besseli(k,x) is the modified Bessel function of the First Kind with order k. 
    // --------------------------------------------------------

    /** 
     * Computes the value of function of PNEBI(0,x) = 2.0 * exp(-x) * besseli(0,x) 
     * using the polynomial approximation of Abramowitz-Stegun (9.8.1)-(9.8.2). 
     * Common library functions computing besseli(0,x) leads to numeric overflow 
     * or exploits inaccurate (and substantially flat in our interval!) Hankel 
     * approximation. 
     */
    double evaluatePnebi0Polynom(double x);

    /** 
     * Evaluates PNEBI function in point x for different orders from 0 to n. 
     * This implementation is based on downward recurring formula as suggested in
     *  
     * Aa Vv, Numerical Recipes in C. The Art of Scientific Computing, 2nd edition, 1992. 
     * 
     * It uses down
     */
    void evaluatePnebiVector(int n, double x, std::vector<double>& pnebis);

    /** 
     * Class for storing a Look-Up Table (LUT) of PNEBI function. 
     * The LUT has adaptive step to achieve a trade-off between memory efficiency 
     * and accuracy. 
     * 
     * The function is defined as:
     * 
     *   PNEBI(k,x) = 2.0 * exp(-x) * besseli(k,x)
     * 
     * where besseli(k,x) is the modified bessel function  of the first kind of 
     * order k computed in x. 
     */
    class PnebiLUT {
    private:

        /** Contains the values of PNEBI functions in x for different orders. 
         */
        struct PnebiPoint {
            double x;
            std::vector<double> y;

            PnebiPoint(int n, double xval) : x(xval), y(n + 1) {
                evaluatePnebiVector(n, x, y);
            }

            bool operator<(const PnebiPoint& pp) const {
                return (x < pp.x);
            }

            bool operator>(const PnebiPoint& pp) const {
                return (x > pp.x);
            }
        };

    public:
        /** 
         * Creates an empty LUT. 
         */
        PnebiLUT();

        /** 
         * Creates PNEBI LUT for order k = 0, .., n and tollerance on the estimated 
         * value of PNEBI function. 
         */
        PnebiLUT(int n, double tol);

        /** 
         * Initializes a new LUT. Warning: it resets previous LUT, if any. 
         */
        void init(int n, double tol);

        /** 
         * Returns the maximum order of PNEBI stored in the LUT.
         */
        int getOrderMax() const {
            return orderMax_;
        }

        /** 
         * Computes the value of PNEBI.
         */
        double eval(int k, double x) const;

        /** 
         * Computes the value of PNEBI for all the available orders.
         */
        void eval(double x, std::vector<double>& y) const;
        
        /**
         * Prints the given number of LUT entries.
         * @param out output stream 
         * @param n number of entries to print (as rows)
         * @param k maximum order of PNEBI to be printed (as columns)
         */
        void printLUT(std::ostream& out,int n,int k = 0);

    private:
        std::vector<PnebiPoint> lut_;
        int orderMax_;
        double tolerance_;
    };
    
    // --------------------------------------------------------
    // INTERVAL FUNCTIONS
    // --------------------------------------------------------

    /** Computes lower and upper bounds of cosine function on a given interval.
     */
    void findLUCos(double a, double b, double& cmin, double& cmax);

    /** Computes lower and upper bounds of Fourier Series (represented by its coefficients)
     * on a given interval.
     * The vector of coefficients coeffs[i] are used in Fourier series:
     *   S(x) = \sum_{i=0}^{n} ( coeffs[2*i] * cos(2*i*x) + coeffs[2*i+1] * sin(2*i*x) )
     */
    void findLUFourier(const std::vector<double>& coeffs, double theta0, double theta1, double& fourierMin, double& fourierfMax);
    
    void fft(const std::vector<double>& funIn, std::vector<double>& coeff, int fourierOrder);

} // end of namespace

#endif //ARS_FUNCTIONS_H