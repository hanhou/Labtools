/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#include "prior.h"

#include <iostream>

void GaussPrior::shrink ( double xmin, double xmax ) {
	double s ( 0.5*(xmax-xmin) ), m ( 0.5*(xmin+xmax) );
	if ( s<std() ) {
		mu = m;
		sg = s;
		var = sg*sg;
		twovar = 2*var;
		rng = GaussRandom ( mu, sg );
		normalization = 1./(sqrt(2*M_PI)*sg);
	}
}

double GaussPrior::ppf ( double p, double start ) const {
	if ( p<=0 || p>=1 )
		throw BadArgumentError ( "Requested probability is outside the range" );
	unsigned int i;
	double x,d;

	if (start==NULL)
		x = mu;
	else
		x = start;

	for ( i=0; i<20; i++ ) {
		d = (cdf ( x ) - p)/pdf ( x );
		x -= d;
		if ( fabs(d) < 1e-7 )
			break;
	}
	return x;
}

void BetaPrior::shrink ( double xmin, double xmax ) {
	double s ( 0.5*(xmax-xmin) ), m ( 0.5*(xmin+xmax) );
	if ( s<std() ) {
		beta = m*(1-m)*(1-m)/(s*s) - 1 + m;
		alpha = m*beta/(1-m);
		normalization = betaf(alpha,beta);
		rng = BetaRandom ( alpha, beta );
	}
}

double BetaPrior::ppf ( double p, double start ) const {
	if ( p<=0 || p>=1 )
		throw BadArgumentError ( "Requested probability is outside the range" );
	double x, d, sx;
	unsigned int i;

	if ( start == NULL )
		x = log ( mean() / ( 1-mean() ) );
	else if ( start<=0 || start>=1 ) {
		throw BadArgumentError ( "Beta Distribution can not be evaluated outside the unit interval" );
	} else
		x = log ( start / (1-start) );

	for ( i=0; i<20; i++ ) {
		sx = 1./(1+exp(-x));
		d = (cdf ( sx ) - p) / ( pdf ( sx ) * sx * (1-sx) );
		x -= d;
		if ( fabs ( d ) < 1e-7 )
			break;
	}

	return 1./(1+exp(-x));
}

void GammaPrior::shrink ( double xmin, double xmax ) {
	double xr ( xmin/xmax );
	k = ((1+xr)/(1-xr));
	k *= k;
	theta = xmax / (k+sqrt(k));
	normalization = pow(theta,k)*exp(gammaln(k));
	rng = GammaRandom ( k, theta );
}

double GammaPrior::ppf ( double p, double start ) const {
	if ( p<=0 || p>=1 )
		throw BadArgumentError ( "Requested probability is outside the range" );
	double x,d;
	unsigned int i;

	if ( start == NULL )
		x = sqrt ( mean() );
	else
		x = sqrt ( start );

	for ( i=0; i<20; i++ ) {
		d = (cdf ( x*x ) - p)/(2*pdf( x*x )*x);
		x -= d;
		if ( fabs ( d ) < 1e-7 )
			break;
	}
	return x*x;
}

void nGammaPrior::shrink ( double xmin, double xmax ) {
	double ymin(-xmax), ymax(-xmin);
	GammaPrior::shrink ( ymin, ymax );
}
