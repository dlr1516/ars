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
#include <ars/functions.h>
#include <cassert> 
#include <algorithm>
#include <unsupported/Eigen/FFT>

#include <ars/definitions.h>


namespace ars {

    // --------------------------------------------------------
    // COS-SIN FAST
    // --------------------------------------------------------

    void fastCosSin(double x, double& c, double& s) {
        constexpr double factor = 1.0 / (2.0 * M_PI);
        x *= factor;

        c = x - (0.25 + floor(x + 0.25));
        c *= 16.0 * (std::abs(c) - 0.5);
        c += 0.225 * c * (std::abs(c) - 1.0);

        s = x - floor(x + 0.5);
        s *= 16.0 * (0.5 - std::abs(s));
        s += 0.225 * s * (std::abs(s) - 1.0);
    }

    double fastAtan(double x) {
        const double a1 = 0.9998660;
        const double a3 = -0.3302995;
        const double a5 = 0.1801410;
        const double a7 = -0.0851330;
        const double a9 = 0.0208351;
        double x2 = x * x;
        double x4 = x2 * x2;
        return x * (a1 + x2 * (a3 + a7 * x4) + x4 * (a5 + a9 * x4));
    }

    double fastAtan2(double x, double y) {
        double ay = fabs(y);
        double ax = fabs(x);
        bool invert = ay > ax;
        double z = invert ? ax / ay : ay / ax; // z in range [0,1]
        double th = fastAtan(z); // th in range [0,M_PI/4] 
        if (invert) th = M_PI_2 - th; // th in range [0,M_PI/2]
        if (x < 0) th = M_PI - th; // th in range [0,M_PI]
        th = copysign(th, y); // th in range [-M_PI,M_PI]
    }

    double evaluateFourier(const std::vector<double>& coeffs, double theta) {
        double val, cth2, sth2, cth, sth, ctmp, stmp;
        int n;

        if (coeffs.size() % 2 != 0) {
            std::cerr << __FILE__ << "," << __LINE__ << ": the number of coefficients must be even: found " << coeffs.size() << std::endl;
        }
        n = (coeffs.size() / 2) - 1;

        cth2 = cos(theta);
        sth2 = sin(theta);
        cth = 1.0;
        sth = 0.0;
        val = 0.0;

        for (int k = 0; k <= n; ++k) {
            val += coeffs[2 * k] * cth + coeffs[2 * k + 1] * sth;
            ctmp = cth2 * cth - sth2 * sth;
            stmp = sth2 * cth + cth2 * sth;
            cth = ctmp;
            sth = stmp;
        }
        return val;
    }

    // --------------------------------------------------------
    // PNEBI FUNCTIONS
    // --------------------------------------------------------

    double evaluatePnebi0Polynom(double x) {
        double t, t2, tinv, val;

        if (x < 0.0) x = -x;
        t = x / 3.75;

        if (t < 1.0) {
            t2 = t*t;
            val = 1.0 + t2 * (3.5156229 + t2 * (3.0899424 + t2 * (1.2067492 + t2 * (0.2659732 + t2 * (0.360768e-1 + t2 * 0.45813e-2)))));
            val = 2.0 * exp(-x) * val;
        } else {
            tinv = 1 / t;
            val = (0.39894228 + tinv * (0.1328592e-1 + tinv * (0.225319e-2 + tinv * (-0.157565e-2 + tinv *
                    (0.916281e-2 + tinv * (-0.2057706e-1 + tinv * (0.2635537e-1 + tinv * (-0.1647633e-1 + tinv * 0.392377e-2))))))));
            val = 2.0 * val / sqrt(x);
        }

        return val;
    }

    void evaluatePnebiVector(int n, double x, std::vector<double>& pnebis) {
        double factor, seqPrev, seqCurr, seqNext;
        if (pnebis.size() < n + 1) {
            pnebis.resize(n + 1);
        }

        if (x < 0.0) x = -x;

        // If x~=0, then BesselI(0,x) = 1.0 and BesselI(k,x) = 0.0 for k > 0.
        // Thus, PNEBI(0,x) = 2.0 and PNEBI(k,x) = 0.0 for k > 0.
        if (x < 1e-6) {
            std::fill(pnebis.begin(), pnebis.end(), 0.0);
            pnebis[0] = 2.0;
            return;
        }

        // Computes bessel function using back recursion
        factor = 2.0 / x;
        seqPrev = 0.0; // bip
        seqCurr = 1.0; // bi
        seqNext = 0.0; // bim
        for (int k = 2 * (n + (int) sqrt(40.0 * n)); k >= 0; --k) {
            seqNext = seqPrev + factor * k * seqCurr;
            seqPrev = seqCurr;
            seqCurr = seqNext;
            if (k <= n) {
                pnebis[k] = seqPrev;
            }
            // To avoid overflow!
            if (seqCurr > BIG_NUM) {
                seqPrev *= SMALL_NUM;
                seqCurr *= SMALL_NUM;
                for (int i = 0; i < pnebis.size(); ++i) {
                    pnebis[i] *= SMALL_NUM;
                }
                //std::cerr << __FILE__ << "," << __LINE__ << ": ANTI-OVERFLOW!" << std::endl;
            }
        }

        double scaleFactor = evaluatePnebi0Polynom(x) / pnebis[0];
        for (int i = 0; i < pnebis.size(); ++i) {
            pnebis[i] = scaleFactor * pnebis[i];
        }
    }

    // --------------------------------------------------------
    // PNEBI LUT
    // --------------------------------------------------------

    PnebiLUT::PnebiLUT() : lut_(), orderMax_(0), tolerance_(0.01) {
        // LUT NOT INITIALIZED!
    }

    PnebiLUT::PnebiLUT(int n, double tol) : lut_(), orderMax_(n), tolerance_(tol) {
        init(n, tol);
    }

    void PnebiLUT::init(int n, double tol) {
        // The LUT computation takes into account the following facts. Let f(k,x) = PNEBI(k,x) (for short).
        // a. f(k,x) <= f(0,x) for all x;
        // b. d/dx I(0,x) = I(1,x) 
        //    -> d/dx f(0,x) = 2.0 d/dx(exp(-x) * I(0,x)) = 2.0 * (-exp(-x) * I(0,x) + exp(-x) * I(1,x))
        //        = 2.0 * exp(-x) * (I(1,x) - I(0,x)) = f(1,x) - f(0,x) = f'(0,x0)
        // c. Locally around point x0, f(0,x) ~= f(0,x0) + (f(1,x0) - f(0,x0)) * (x - x0)
        //    
        // Property a) grants that, if f(0,x) is sampled according to the desired accuracy, all f(k,x)
        // are sampled with at least the same accuracy. 
        //
        // Function f(0,x) is strictly decreasing w.r.t. x. Thus, f'(0,x) < 0.0 for all x. 
        // Moreover, f'(0,x) is bounded from -2.0 to 0.0 (for x -> +inf). 
        // Suppose to sample f(0,x) in x_i (x_i < x_{i+1} for all i) s.t.  
        //    f(0,x_i) - f(0,x_{i+1}) = toll > 0
        // Then, the function is approximated by a line:
        //    f(0,x_{i+1}) = f(0,x_i) + f'(0,x_i) * (x_{i+1} - x_i) + r(x_{i+1}-x_i)  
        //    x_{i+1} - x_i = (f(0,x_{i+1}) - f(0,x_i)) / f'(0,x_i) - r(x_{i+1}-x_i) 
        //    x_{i+1} = x_i - toll / f'(0,x_i) - r(x_{i+1}-x_i)
        // (the term r(x_{i+1}-x_i) >= 0 and is O((x_{i+1}-x_i)^2)).
        //
        // Case  f'(0,x0) = 0 occurs only for x0 -> +inf. When |f'(0,x)| < eps, the algorithm stops.
        int lutNum = (int) ceil(2.0 / tol) + 20;
        double x, y0prev, ratio, deriv;

        // Inserts initial point x=0.0
        x = 0.0;
        orderMax_ = n;
        tolerance_ = tol;
        lut_.clear();
        lut_.reserve(lutNum);
        lut_.push_back(PnebiPoint(n, x));
        deriv = lut_.back().y[1] - lut_.back().y[0]; // Initial value of deriv is -2.0
        y0prev = lut_.back().y[0];
        //  std::cout << "x " << x << ", deriv " << deriv << std::endl;
        for (int i = 1; i < lutNum; ++i) {
            x = x - tol / deriv;
            lut_.push_back(PnebiPoint(n, x));
            ratio = deriv * (y0prev - lut_.back().y[0]) / tol;
            deriv = std::min(std::max(lut_.back().y[1] - lut_.back().y[0], ratio), -10.0 * std::numeric_limits<double>::min());
            //            std::cout << "x " << x << ", y[0] " << lut_.back().y[0] << " (incr prev y[0] " << (lut_.back().y[0] - y0prev) << "), "
            //                      << "real deriv " << (lut_.back().y[1] - lut_.back().y[0]) << ", deriv " << deriv << ", ratio " << ratio << std::endl;
            y0prev = lut_.back().y[0];
        }
        //std::cout << __FILE__ << "," << __LINE__ << ": lut num " << lut_.size() << ", x max " << lut_.back().x << " -> y[0] " << lut_.back().y[0] << std::endl;
    }

    double PnebiLUT::eval(int k, double x) const {
        double val, w;

        if (lut_.empty()) {
            std::cerr << __FILE__ << "," << __LINE__ << ": empty LUT" << std::endl;
            return 0.0;
        }

        // Checks that 1) argument x is positive; 2) order of PNEBI has been computed;
        if (x < 0.0) x = -x;
        assert(0 <= k && k < orderMax_);

        PnebiPoint tmp(0, x);
        auto upper = std::upper_bound(lut_.begin(), lut_.end(), tmp);
        auto lower = upper;
        std::advance(lower, -1);

        if (upper == lut_.end()) {
            //std::cout << __FILE__ << "," << __LINE__ << ": lower->x " << lower->x << std::endl;
            val = lower->y[k];
        } else {
            //    std::cout << __FILE__ << "," << __LINE__ << ": lower->x " << lower->x << ", upper->x " << upper->x << std::endl;
            w = (x - lower->x) / (upper->x - lower->x);
            val = (1.0 - w) * lower->y[k] + w * upper->y[k];
        }
        return val;
    }

    void PnebiLUT::eval(double x, std::vector<double>& y) const {
        double val, w;

        // Checkes LUT initialization
        if (lut_.empty()) {
            std::cerr << __FILE__ << "," << __LINE__ << ": empty LUT" << std::endl;
            return;
        }

        // Checks that 1) argument x is positive; 2) order of PNEBI has been computed;
        if (x < 0.0) x = -x;

        // Checkes the size of vector (and adapt it if needed)
        if (y.size() < orderMax_ + 1) {
            y.resize(orderMax_ + 1);
        }

        PnebiPoint tmp(0, x);
        auto upper = std::upper_bound(lut_.begin(), lut_.end(), tmp);
        auto lower = upper;
        std::advance(lower, -1);

        if (upper == lut_.end()) {
            std::copy(lower->y.begin(), lower->y.end(), y.begin());
            //evaluatePnebiVector(orderMax_,x,y);
        } else {
            w = (x - lower->x) / (upper->x - lower->x);
            for (int i = 0; i < lower->y.size() && i < upper->y.size() && i < y.size(); ++i) {
                y[i] = (1.0 - w) * lower->y[i] + w * upper->y[i];
            }
        }
    }

    void PnebiLUT::printLUT(std::ostream& out, int n, int k) {
        int incr = std::max<int>(lut_.size() / n, 1);
        for (int i = 0; i < lut_.size(); i += incr) {
            out << i << "\t" << lut_[i].x << "\t";
            for (int j = 0; j <= k && j < lut_[i].y.size(); ++j) {
                out << lut_[i].y[j] << "\t";
            }
            out << std::endl;
        }
        out << (lut_.size() - 1) << "\t" << lut_.back().x << "\t";
        for (int j = 0; j <= k && j < lut_.back().y.size(); ++j) {
            out << lut_.back().y[j] << "\t";
        }
        out << std::endl;
    }

    // --------------------------------------------------------
    // INTERVAL FUNCTIONS
    // --------------------------------------------------------

    void findLUCos(double a, double b, double& cmin, double& cmax) {
        double amod, bmod;

        if (a > b) std::swap(a, b);

        if (b - a >= 2.0 * M_PI) {
            cmin = -1.0;
            cmax = +1.0;
        } else {
            // Normalizes to circular interval [0, 2*M_PI[
            amod = fmod(a, 2.0 * M_PI);
            if (amod < 0.0) amod += 2.0 * M_PI;
            bmod = fmod(b, 2.0 * M_PI);
            if (bmod < 0.0) bmod += 2.0 * M_PI;
            // Case bmod < amod: for example [300,30[ deg: angle 0 is included.
            if (bmod < amod) {
                cmax = +1.0;
                if (bmod < M_PI && M_PI < amod) {
                    cmin = std::min(cos(amod), cos(bmod));
                } else {
                    cmin = -1.0;
                }
                //      if (M_PI < bmod || amod < M_PI) {
                //        cmin = -1.0;
                //      }
                //      else {
                //        cmin = std::min(cos(amod),cos(bmod));
                //      }
            } else {
                cmax = std::max(cos(amod), cos(bmod));
                if (amod < M_PI && M_PI < bmod) {
                    cmin = -1.0;
                } else {
                    cmin = std::min(cos(amod), cos(bmod));
                }
            }
        }
    }

    void findLUFourier(const std::vector<double>& coeffs, double theta0, double theta1, double& fourierMin, double& fourierMax) {
        double amplitude, phase, sinusoidMin, sinusoidMax;
        int n, i0, i1;

        if (coeffs.size() % 2 != 0) {
            std::cerr << __FILE__ << "," << __LINE__ << ": the number of coefficients must be even: found " << coeffs.size() << std::endl;
        }
        n = (coeffs.size() / 2) - 1;

        if (theta1 < theta0) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid interval [" << theta0 << "," << theta1 << "]: swapping endpoints to continue" << std::endl;
            std::swap(theta0, theta1);
        }

        // fourierMin and fourierMax initialized with constant component
        fourierMin = coeffs[0];
        fourierMax = coeffs[0];
        for (int k = 1; k <= n; ++k) {
            // t_k = a_k * cos(2*k*theta) + b_k * sin(2*k*theta) = amplitude * cos(2*k*theta - phase)
            // Period of k-th terms is M_PI / k. 
            amplitude = sqrt(coeffs[2 * k] * coeffs[2 * k] + coeffs[2 * k + 1] * coeffs[2 * k + 1]);
            phase = atan2(coeffs[2 * k + 1], coeffs[2 * k]);
            //std::cout << "k " << k << ", amplitude " << amplitude << ", phase[deg] " << (180.0/M_PI*phase) << std::endl;
            // If the [theta0,theta1] interval is larger than period, then the whole sinusoid amplitude is considered.
            // Otherwise, a more refined evaluation is performed.
            findLUCos(2 * k * theta0 - phase, 2 * k * theta1 - phase, sinusoidMin, sinusoidMax);
            fourierMin += amplitude * sinusoidMin;
            fourierMax += amplitude * sinusoidMax;
        }
    }

    void fft(const std::vector<double>& funIn, std::vector<double>& coeffs, int fourierOrder) {
        Eigen::FFT<double> fft_;
        std::vector<std::complex<double> > freqvec;
        int n = funIn.size();

        fft_.fwd(freqvec, funIn);
        
        ARS_PRINT("fft: input size " << funIn.size() << ", output complex size " << freqvec.size());
        
        coeffs.resize(2 * n + 2);
        coeffs[0] = freqvec[0].real();
        coeffs[1] = 0.0;
        for (int i = 1; i < n && i < fourierOrder; ++i) {
            coeffs[2 * i] = 0.5 * (freqvec[i].real() + freqvec[n - i].real());
            coeffs[2 * i + 1] = 0.5 * (freqvec[i].imag() - freqvec[n - i].imag());
        }
    }

} // end of namespace

