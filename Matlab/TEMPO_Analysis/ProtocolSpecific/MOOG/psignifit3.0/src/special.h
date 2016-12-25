/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#ifndef SPECIAL_H
#define SPECIAL_H

#include <cmath>
#include <cstdlib>

/** \brief gaussian cumulative distribution function */
double Phi ( double x );

/** \brief inverse of the gaussian cumulative distribution function
 *
 * This function is really expensive. It inverts the gaussian cdf by performing a numerical
 * solution of the equation
 * Phi(x) - p = 0
 */
double invPhi ( double p );

/** \brief logarithm that does not return nan but instead a very low value (-1e20) */
double safe_log ( double x );

/** \brief logarithm of the gamma function */
double gammaln ( double xx );

/** \brief incomplete unnormalized gamma function */
double gammainc ( double x, double a );

/** beta function */
double betaf ( double z, double w );

/** regularized incomplete beta function based on continued fractions */
double betainc ( double x, double al, double bt );

/** psi function, i.e. d log Gamma / dx */
double psi ( double z );

/** digamma (derivative of psi function) */
double digamma ( double z );

#endif
