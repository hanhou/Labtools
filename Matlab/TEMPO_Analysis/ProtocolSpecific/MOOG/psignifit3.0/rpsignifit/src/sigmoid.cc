/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#include "sigmoid.h"

#ifdef DEBUG_SIGMOID
#include <iostream>
#endif

/** Logistic Sigmoid **********************************************************/

double PsiLogistic::f ( double x ) const
{
#ifdef DEBUG_SIGMOID
	std::cerr << "In logistic function.\n";
#endif
	/*
	if (x!=lastx) {
		lastx = x;
		lastfx = 1./(1.+exp(-x));
	}
	*/
#ifdef DEBUG_SIGMOID
	std::cerr << " lastx = " << lastx << "\n lastfx = " << lastfx << "\n";
#endif
	return 1./(1.+exp(-x));
}

double PsiLogistic::df ( double x ) const
{
	return f(x)*(1-f(x));
}

double PsiLogistic::ddf ( double x ) const
{
	return f(x)*(1-f(x))*(1-2*f(x));
}

/** Gauss Sigmoid *************************************************************/

double PsiGauss::f ( double x ) const
{
	/*
	if (x==lastx)
		return lastf;
	else {
		lastx = x;
		lastf = Phi(x);
		return lastf;
	}
	*/
	return Phi(x);
}

double PsiGauss::df ( double x ) const
{
	/*
	if (x==lastx_d)
		return lastdf;
	else {
		lastx_d = x;
		lastdf = exp ( - 0.5*x*x ) / sqrt(2*M_PI);
		return lastdf;
	}
	*/
	return exp ( - 0.5*x*x ) / sqrt(2*M_PI);
}

double PsiGauss::ddf ( double x ) const
{
	/*
	if (x==lastx_dd)
		return lastddf;
	else {
		lastx_dd = x;
		lastddf = -x*df(x);
		return lastddf;
	}
	*/
	return -x*df(x);
}

double PsiGauss::inv ( double p ) const
{
	/*
	if (p==lastp)
		return lastinvp;
	else {
		lastp = p;
		lastinvp = invPhi(p);
		return lastinvp;
	}
	*/
	return invPhi(p);
}

/** Gumbel_l Sigmoid *********************************************************/

double PsiGumbelL::f ( double x ) const
{
	/*
	if (x!=lastx) {
		lastx = x;
		lastf = 1-exp(-exp(x));
	}
	return lastf;
	*/
	return 1-exp(-exp(x));
}

double PsiGumbelL::df ( double x ) const
{
	/*
	if ( x!=lastdx ) {
		lastdx = x;
		lastdf = exp( x - exp(x));
	}
	return lastdf;
	*/
	return exp( x - exp(x));
}

double PsiGumbelL::ddf ( double x ) const
{
	/*
	if ( x!=lastddx ) {
		lastddx = x;
		lastddf = exp ( x - exp(x) ) * (1-exp(x));
	}
	return lastddf;
	*/
	return exp ( x - exp(x) ) * (1-exp(x));
}

double PsiGumbelL::inv ( double p ) const
{
	/*
	if ( p!=lastp ) {
		lastp = p;
		lastinvp = log(-log(1-p));
	}
	return lastinvp;
	*/
	return log(-log(1-p));
}

/** Gumbel_r Sigmoid **********************************************************/

double PsiGumbelR::f ( double x ) const
{
	/*
	if (x!=lastx) {
		lastx = x;
		lastf = exp(-exp(-x));
	}
	return lastf;
	*/
	return exp(-exp(-x));
}

double PsiGumbelR::df ( double x ) const
{
	/*
	if ( x!=lastdx ) {
		lastdx = x;
		lastdf = exp(-x-exp(-x));
	}
	return lastdf;
	*/
	return exp(-x-exp(-x));
}

double PsiGumbelR::ddf ( double x ) const
{
	/*
	if ( x!=lastddx ) {
		lastddx = x;
		lastddf = exp ( -x - exp(-x) ) * (exp(-x)-1);
	}
	return lastddf;
	*/
	return exp ( -x - exp(-x) ) * (exp(-x)-1);
}

double PsiGumbelR::inv ( double p ) const
{
	/*
	if ( p!=lastp ) {
		lastp = p;
		lastinvp = -log(-log(p));
	}
	return lastinvp;
	*/
	return -log(-log(p));
}

/** Cauchy Sigmoid ************************************************************/

double PsiCauchy::f ( double x ) const
{
	return atan ( x )/M_PI + 0.5;
}

double PsiCauchy::df ( double x ) const
{
	return 1./(M_PI*(1+x*x));
}

double PsiCauchy::ddf ( double x ) const
{
	return -2*x/( M_PI * (1+2*x*x+x*x*x*x) );
}

double PsiCauchy::inv ( double p ) const
{
	return tan ( M_PI*(p-0.5) );
}

/** Exponential cdf ***********************************************************/

double PsiExponential::f ( double x ) const
{
	if (x<0)
		return 0;
	else
		return 1-exp ( -x );
}

double PsiExponential::df ( double x ) const
{
	if (x<0)
		return 0;
	else
		return exp( -x );
}

double PsiExponential::ddf ( double x ) const
{
	if (x<0)
		return 0;
	else
		return -exp( -x );
}

double PsiExponential::inv ( double p ) const throw(BadArgumentError)
{
	if ( p>0 && p<1 )
		return -log(1-p);
	else
		throw BadArgumentError("PsiExponential.inv is only valid in the range 0<x<1");
}

