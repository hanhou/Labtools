/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#ifndef RNG_H
#define RNG_H

#include <cstdlib>
#include <cmath>
#include "errors.h"


class PsiRandom
{
	public:
		PsiRandom ( void ) {}
		double rngcall ( void );
		virtual double draw ( void ) { throw NotImplementedError(); }
		virtual PsiRandom * clone ( void ) const {throw NotImplementedError(); }
};

class GaussRandom : public PsiRandom
{
	private:
		double mu;
		double sigma;
		bool good;
		double x1;
		double x2;
		double w;
		double y;
	public:
		GaussRandom ( double mean=0, double standarddeviation=1 ) : mu ( mean ), sigma ( standarddeviation ), good ( false ) {}
		double draw ( void );              ///< draw a random number using box muller transform
		PsiRandom * clone ( void ) const { return new GaussRandom(*this); }
};

class UniformRandom : public PsiRandom
{
	private:
		double lower;
		double upper;
	public:
		UniformRandom ( double low=0, double up=1 ) : lower(low), upper(up) {}
		double draw ( void ) { return  (upper-lower)*rngcall()+lower; }
		PsiRandom * clone ( void ) const { return new UniformRandom(*this); }
};

class BinomialRandom : public PsiRandom
{
	private:
		int n;
		double p;
	public:
		BinomialRandom ( int number, double probability ) : n(number), p(probability) {}
		double draw ( void );
		void setprm ( int number, double probability ) { n = number; p = probability; }
		PsiRandom * clone ( void ) const { return new BinomialRandom(*this); }
};

class GammaRandom : public PsiRandom
{
	private:
		double k;
		double theta;
                GaussRandom grng;
	public:
		GammaRandom ( double shape, double scale ) : k (shape), theta(scale), grng() {}
		double draw ( void );              ///< draw a random number
		PsiRandom * clone ( void ) const { return new GammaRandom(*this); }
};

class BetaRandom : public PsiRandom
{
	private:
		double alpha;
		double beta;
                GammaRandom grnga;
                GammaRandom grngb;
	public:
		BetaRandom ( double alpha, double beta ) : alpha(alpha), beta(beta), grnga (alpha, 1), grngb (beta, 1) {}
		double draw ( void );              ///< draw a random number
		PsiRandom * clone ( void ) const { return new BetaRandom(*this); }
};


void setSeed(long int seedval);

#endif
