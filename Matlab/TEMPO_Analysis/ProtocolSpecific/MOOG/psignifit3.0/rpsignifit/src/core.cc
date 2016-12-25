/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#include "core.h"

/************************************************************
 * abCore methods
 */

double abCore::dg ( double x, const std::vector<double>& prm, int i ) const {
	switch (i) {
	case 0:
		return -1./prm[1];
		break;
	case 1:
		return -(x-prm[0])/(prm[1]*prm[1]);
		break;
	default:
		// If the parameter does not exist in the abCore the derivative with respect to it will always be 0
		return 0;
		break;
	}
}

double abCore::dgx ( double x, const std::vector<double>& prm ) const {
	return 1./prm[1];
}

double abCore::ddg ( double x, const std::vector<double>& prm, int i, int j ) const {
	if (i==j) {
		switch (i) {
		case 0:
			return 0;
			break;
		case 1:
			return 2*(x-prm[0])/(prm[1]*prm[1]*prm[1]);
			break;
		default:
			// If the parameter does not exist in the abCore the derivative with respect to it will always be 0
			return 0;
			break;
		}
	} else if ((i==0 && j==1) || (i==1 && j==0)) {
	    return 1./(prm[1]*prm[1]);
	} else
		    // If the parameter does not exist in the abCore the derivative with respect to it will always be 0
		return 0;
}

double abCore::inv ( double y, const std::vector<double>& prm ) const {
	return y*prm[1] + prm[0];
}

double abCore::dinv ( double y, const std::vector<double>& prm, int i ) const {
	switch (i) {
	case 0:
		return 1;
		break;
	case 1:
		return y;
		break;
	default:
		return 0;
		break;
	}
}

std::vector<double> abCore::transform ( int nprm, double a, double b ) const {
	std::vector<double> out ( nprm, 0 );
	out[1] = 1./b;
	out[0] = -a/b;
	return out;
}

/************************************************************
 * mwCore methods
 */
mwCore::mwCore( const PsiData* data, const int sigmoid, const double alpha )
	: sigmtype(sigmoid), alpha(alpha), zshift(0) {
	switch (sigmoid) {
		case 1:
			// logistic
			zalpha = 2*log(1./alpha-1.);
			break;
		case 2:
			// gaussian
			zalpha = invPhi(1-alpha)-invPhi(alpha);
			break;
		case 3:
			// gumbel
			zalpha = log(-log(alpha))-log(-log(1.-alpha));
			zshift = log(-log(0.5));
			break;
		case 4:
			// cauchy
			zalpha = -2*tan(M_PI*(alpha-0.5));
			zshift = 0;
			break;
		case 5:
			// Exponential
			zalpha = log( (1-alpha)/alpha );
			zshift = log(2.);
			break;
		case 6:
			// gumbel_r
			zalpha = -log(-log(1.-alpha))+log(-log(alpha));
			zshift = -log(-log(0.5));
			break;
		default:
			throw NotImplementedError();
	}
}

double mwCore::g ( double x, const std::vector<double>& prm ) const {
	return zalpha*(x-prm[0])/prm[1] + zshift;
}

double mwCore::dg ( double x, const std::vector<double>& prm, int i ) const {
	switch (i) {
		case 0:
			return -zalpha/prm[1];
			break;
		case 1:
			return -zalpha*(x-prm[0])/(prm[1]*prm[1]);
			break;
		default:
			// Irrelevant parameter ~> derivative is 0.
			return 0;
			break;
	}
}

double mwCore::dgx ( double x, const std::vector<double>& prm ) const {
	return zalpha/prm[1];
}

double mwCore::ddg ( double x, const std::vector<double>& prm, int i, int j ) const {
	if (i==j) {
		if (i==0)
			return 0;
		else if (i==1)
			return 2*zalpha*(x-prm[0])/(prm[1]*prm[1]*prm[1]);
		else
			return 0;
	} else if ( (i==0 && j==1) || (i==1 && j==0) )
		return zalpha/(prm[1]*prm[1]);
	else
		return 0;
}

double mwCore::inv ( double y, const std::vector<double>& prm ) const {
	return prm[0] + prm[1]*(y-zshift)/zalpha;
}

double mwCore::dinv ( double p, const std::vector<double>& prm, int i ) const {
	switch (i) {
		case 0:
			return 1;
			break;
		case 1:
			return (p-zshift)/zalpha;
			break;
		default:
			return 0;
			break;
	}
}

std::vector<double> mwCore::transform ( int nprm, double a, double b ) const {
	std::vector<double> out ( nprm, 0 );
	out[1] = zalpha/b;
	out[0] = out[1]*(zshift-a)/zalpha;
	return out;
}

/************************************************************
 * logarithmicCore
 */

double logCore::g ( double x, const std::vector<double>& prm ) const throw(BadArgumentError)
{
	if (x<0)
		throw BadArgumentError("logCore.g is only valid in the range x>=0");
	return prm[0] * (x==0 ? -1e10 : log(x)) + prm[1];
}

logCore::logCore( const PsiData* data, const int sigmoid, const double alpha ) : scale(0) {
	unsigned int i;
	// we need this to scale starting values obtained from logistic regression so that they are correct "on average"
	for (i=0; i<data->getNblocks(); i++)
		scale += data->getIntensity(i)/log(data->getIntensity(i));
	scale /= data->getNblocks();
}

double logCore::dg ( double x, const std::vector<double>& prm, int i ) const {
	switch (i) {
		case 0:
			return log(x);
			break;
		case 1:
			return 1;
			break;
		default:
			return 0;
			break;
	}
}

double logCore::dgx ( double x, const std::vector<double>& prm ) const {
	return prm[0]/x;
}


double logCore::dinv ( double y, const std::vector<double>& prm, int i ) const {
	switch (i) {
		case 0:
			return exp((y-prm[1])/prm[0]) * (prm[1]-y)/(prm[0]*prm[0]);
			break;
		case 1:
			return -exp((y-prm[1])/prm[0])/prm[0];
			break;
		default:
			return 0;
			break;
	}
}

std::vector<double> logCore::transform ( int nprm, double a, double b ) const {
	std::vector<double> prm ( nprm, 0 );
	prm[0] = b*scale;  // we scale the intercept so that it is correct "on average"
	prm[1] = a;
	return prm;
}

/************************************************************
 * weibullCore
 */

weibullCore::weibullCore( const PsiData* data, const int sigmoid, const double alpha ) : twooverlog2(2./log(2)) , loglog2 ( log(log(2.)) )
{
	// approximate log(x)~ax+b by a linear function over the range of x values in data
	double covxlogx(0),varx(0);
	double meanx(0), meanlogx(0);
	unsigned int i;
	for (i=0; i<data->getNblocks(); i++) {
		meanx += data->getIntensity(i);
		meanlogx += log(data->getIntensity(i));
	}

	meanx /= data->getNblocks();
	meanlogx /= data->getNblocks();

	for (i=0; i<data->getNblocks(); i++) {
		varx     += pow( data->getIntensity(i)-meanx, 2);
		covxlogx += (data->getIntensity(i)-meanx) * (log(data->getIntensity(i))-meanlogx);
	}
	varx     /= data->getNblocks()-1;
	covxlogx /= data->getNblocks()-1;

	loglina = covxlogx/varx;
	loglinb = meanlogx - loglina*meanx;
}

double weibullCore::dg ( double x, const std::vector<double>& prm, int i ) const throw(BadArgumentError)
{
	if (x<0)
		throw BadArgumentError("weibullCore.dg is only valid in the range x>=0");

	if (i==0) {
		// return -twooverlog2*prm[1] * log(prm[0]);
		return twooverlog2*prm[1] * ( log(x) - log(prm[0]) - 1);
	} else if (i==1) {
		return twooverlog2*prm[0]*((x==0 ? -1e10 : log(x))-log(prm[0]));
	} else {
		return 0;
	}
}

double weibullCore::dgx ( double x, const std::vector<double>& prm ) const {
	return twooverlog2*prm[0]*prm[1]/x;
}

double weibullCore::ddg ( double x, const std::vector<double>& prm, int i, int j ) const throw(BadArgumentError)
{
	if (x<0)
		throw BadArgumentError("weibullCore.ddg is only valid in the range x>=0");

	if (i==j) {
		if (i==0)
			return -twooverlog2 * prm[1] / prm[0];
		else
			return 0;
	} else {
		if ( (i==0 && j==1) || (i==1 && j==0) ) {
			// return -twooverlog2 * log(prm[0]);
			return twooverlog2 * ( log(x) - log(prm[0]) - 1 );
		} else
			return 0;
	}
}

double weibullCore::inv ( double y, const std::vector<double>& prm ) const
{
	return prm[0] * exp (y/(prm[0]*prm[1]*twooverlog2));
}

double weibullCore::dinv ( double y, const std::vector<double>& prm, int i ) const
{
	if (i==0)
		return exp(y/(prm[0]*prm[1]*twooverlog2)) * (1-y/(twooverlog2*prm[0]*prm[1]));
	else if (i==1)
		return -exp(y/(prm[0]*prm[1]*twooverlog2))*y/(twooverlog2*prm[1]*prm[1]);
	else
		return 0;
}

std::vector<double> weibullCore::transform ( int nprm, double a, double b ) const
{
	std::vector<double> prm ( nprm, 0 );
	prm[1] = exp ( b/loglinb );
	prm[0] = ((a/loglinb)/twooverlog2)/prm[1];

	return prm;
}

/************************************************************
 * polyCore
 */

polyCore::polyCore( const PsiData* data, const int sigmoid, const double alpha )
{
	double meanx (0),varx(0);
	unsigned int i;

	for (i=0; i<data->getNblocks(); i++) {
		meanx += data->getIntensity(i);
	}
	meanx /= data->getNblocks();

	for (i=0; i<data->getNblocks(); i++) {
		varx += pow( data->getIntensity(i)-meanx, 2 );
	}
	varx /= data->getNblocks();
	varx = sqrt(varx);

	x1 = meanx+varx;
	x2 = meanx-varx;
}

double polyCore::dg ( double x, const std::vector<double>& prm, int i ) const
{
	if (x<0)
		return 0;
	else {
		if (i==0)
			return -prm[1] * x * pow(x/prm[0] , prm[1]-1)/(prm[0]*prm[0]);
		else if (i==1)
			return pow (x/prm[0], prm[1] ) * log(x/prm[0]);
		else
			return 0;
	}
}

double polyCore::dgx ( double x, const std::vector<double>& prm ) const
{
	if (x<0)
		return 0;
	else {
		return prm[1]*pow(prm[0],-prm[1])*pow(x,prm[1]-1);
	}
}

double polyCore::ddg ( double x, const std::vector<double>& prm, int i, int j ) const
{
	if (x<0)
		return 0;
	else {
		if (i==j) {
			if (i==0)
				return prm[1]*x*(prm[1]+1)*pow(x/prm[0],prm[1]-1)/(prm[0]*prm[0]*prm[0]);
			else if (i==1)
				return pow(x/prm[0],prm[1]) * pow(log(x/prm[0]),2);
			else
				return 0;
		} else if ( (i==0 && j==1) || (j==0 && i==1) ) {
			return - pow(x/prm[0],prm[1]) * ( prm[1]*log(x/prm[0]) + 1. ) / prm[0];
		} else
			return 0;
	}
}

double polyCore::inv ( double y, const std::vector<double>& prm ) const
{
	return prm[0] * pow ( y, 1./prm[1] );
}

double polyCore::dinv ( double y, const std::vector<double>& prm, int i ) const
{
	if (i==0) {
		return pow ( y, 1./prm[1] );
	} else if (i==1) {
		return - log(y) * prm[0] * pow ( y, 1./prm[1] )/(prm[1]*prm[1]);
	} else
		return 0;
}

std::vector<double> polyCore::transform ( int nprm, double a, double b ) const
{
	std::vector<double> prm ( nprm, 0 );

	if ( a+b*x1 < 0 )
		a = -b*x1+.1;
	if ( a+b*x2 < 0 )
		a = -b*x2+.1;

	prm[1] = log ( (a+b*x2)/(a+b*x1) )/log(x2/x1);
	prm[0] = x1*pow(a+b*x1, - 1./prm[1]);

	return prm;
}

/************************************************************
 * NakaRushton
 */

NakaRushton::NakaRushton ( const PsiData *data, const int sigmoid, const double alpha )
	: x(data->getNblocks())
{
	unsigned int i;

	for ( i=0; i<data->getNblocks(); i++ ) {
		x[i] = data->getIntensity(i);
	}
}

double NakaRushton::dg ( double x, const std::vector<double>& prm, int i ) const
{
	double sigm, k;
	double xk,sigmk;
	if (x<0)
		return 0;
	sigm = prm[0];
	k    = prm[1];
	xk = pow(x,k);
	sigmk = pow(sigm,k);

	// These derivatives are from sympy
	switch (i) {
		case 0:
			return -k * xk * sigmk / ( sigm * pow( xk + sigmk, 2) );
			break;
		case 1:
			return xk*log(x)/(xk+sigmk) - xk * ( xk * log(x) + sigmk*log(sigm) ) / pow(xk+sigmk,2);
			break;
		default:
			return 0;
	}
}

double NakaRushton::ddg ( double x, const std::vector<double>& prm, int i, int j ) const
{
	double sigm,k;
	double xk,sigmk;
	double logx,logsigm;
	if (x<0)
		return 0;
	sigm = prm[0];
	k    = prm[1];
	xk   = pow(x,k);
	sigmk= pow(sigm,k);
	logx = log(x);
	logsigm = log(sigm);

	// These derivatives are from sympy
	if ( (i==0) && (j==0) ) {
		return 2*xk*k*k*sigmk*sigmk/(sigm*sigm*pow(xk+sigmk,3))
			+ (k*xk*sigmk - xk*k*k*sigmk)/(sigm*sigm*pow(xk+sigmk,2));
	} else if ( (i==1) && (j==1) ) {
		return -xk * (xk*logx*logx + sigmk*logsigm*logsigm)/pow(xk+sigmk,2)
			+ xk*(xk*logx+sigmk*logsigm)*(2*xk*logx+2*sigmk*logsigm)/pow(xk+sigmk,3)
			- 2*xk*(xk*logx+sigmk*logsigm)*logx/pow(xk+sigmk,2)
			+ xk*logx*logx/(xk+sigmk);
	} else if ( ((i==0) && (j==1)) || ((i==1) && (j==0)) ) {
		return -xk*(x*sigmk*logsigm + xk)/(sigm*pow(xk+sigmk,2))
			- k*xk*sigmk*logx / (sigm*pow(xk+sigmk,2))
			+ 2*k*xk*sigmk*(xk*logx+sigmk*logsigm)/(sigm*pow(xk+sigmk,3));
	} else
		return 0;
}

double NakaRushton::dgx ( double x, const std::vector<double>& prm ) const
{
	double xmk ( pow(x,prm[1]-1) );
	double xk (xmk*x), sgk ( pow(prm[0], prm[1]) );
	if ( x<0 )
		return 0;
	else
		return prm[1]*xmk*(sgk+xk+xk)/((sgk+xk)*(sgk+xk));
}

double NakaRushton::inv ( double y, const std::vector<double>& prm ) const
{
	return pow ( pow(prm[0], prm[1])*y/(1-y), 1.0/prm[1] );
}

double NakaRushton::dinv ( double y, const std::vector<double>& prm, int i ) const
{
	/*
	double sqrtpart, k, sigm, logsigm, sigmk;
	sigmk    = pow(sigm,k);
	sqrtpart = pow ( sigmk / (1-y), 1.0/prm[1] );
	logsigm  = log ( sigm );
	return sqrtpart * ( logsigm/k - log(sigmk/(1-y))/(k*k));
	*/
	double odds (y/(1-y)), sigm(prm[0]), k(prm[1]);
	switch (i) {
		case 0:
			return pow(odds, 1./k);
			break;
		case 1:
			return sigm*pow(odds, 1./k) * (
					log(sigm)/k - log(odds*pow(sigm,k))/(k*k)
					);
	}
	return -1;
}

std::vector<double> NakaRushton::transform ( int nprm, double a, double b ) const
{
	double s1(0),s2(0),s3(0),s4(0), logxi, xi;
	double khat, klogsigmhat;
	unsigned int i;

	// Using linear regression to approximate the logarithm linearly on the desired stimulus range
	for ( i=0; i<x.size(); i++ ) {
		xi    = x[i];
		logxi = log(xi);
		s1 += logxi*(a+b*xi);
		s2 += logxi;
		s3 += a+b*xi;
		s4 += logxi*logxi;
	}

	khat = s1-s2*s3;
	khat /= (s4-s2*s2);
	s2 /= x.size();
	s3 /= x.size();
	klogsigmhat = s3 - khat * s2;

	std::vector<double> prm ( nprm );
	prm[1] = khat;
	prm[0] = exp ( klogsigmhat/khat );

	return prm;
}
