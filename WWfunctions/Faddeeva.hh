/* C++ declarations for the Faddeeva library (Steven G. Johnson, MIT).
   This header declares the Faddeeva::  namespace used by Faddeeva_impl.hh.
   The C API is in Faddeeva.h; do NOT include that from C++ code (it pulls in
   <complex.h> and uses C99 "double complex", both of which break under Cling). */

#ifndef FADDEEVA_HH
#define FADDEEVA_HH 1

#include <complex>

namespace Faddeeva {

using std::complex;

// w(z) = exp(-z^2) erfc(-iz)  [Faddeeva / scaled complex error function]
extern complex<double> w(complex<double> z, double relerr = 0);
extern double w_im(double x);          // Im[w(x)] for real x

// erfcx(z) = exp(z^2) erfc(z)
extern complex<double> erfcx(complex<double> z, double relerr = 0);
extern double erfcx(double x);

// erf(z)
extern complex<double> erf(complex<double> z, double relerr = 0);
extern double erf(double x);

// erfi(z) = -i erf(iz)
extern complex<double> erfi(complex<double> z, double relerr = 0);
extern double erfi(double x);

// erfc(z) = 1 - erf(z)
extern complex<double> erfc(complex<double> z, double relerr = 0);
extern double erfc(double x);

// Dawson(z) = sqrt(pi)/2 * exp(-z^2) * erfi(z)
extern complex<double> Dawson(complex<double> z, double relerr = 0);
extern double Dawson(double x);

} // namespace Faddeeva

#endif // FADDEEVA_HH
