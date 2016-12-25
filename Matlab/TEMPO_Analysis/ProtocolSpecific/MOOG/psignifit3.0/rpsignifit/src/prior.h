/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#ifndef PRIOR_H
#define PRIOR_H

#include <cmath>
#include "rng.h"
#include "special.h"

/** \brief base class for all priors
 *
 * This default prior does nothing in particular. It poses no restriction on the respective parameter at all.
 * at any value x, the prior returns 1. Thus it is an "improper" prior in the sense that it does not correspond
 * to a proper probability distribution.
 */
class PsiPrior
{
	private:
		PsiRandom rng;
	public:
		virtual double pdf ( double x ) const { return 1.;}    ///< evaluate the pdf of the prior at position x (in this default form, the parameter is completely unconstrained)
		virtual double dpdf ( double x ) { return 0.; }  ///< evaluate the derivative of the pdf of the prior at position x (in this default form, the parameter is completely unconstrained)
		virtual double rand ( void ) { return rng.draw(); } ///< draw a random number
		virtual PsiPrior * clone ( void ) const { throw NotImplementedError(); }///< clone by value
		virtual double mean ( void ) const { return 0; } ///< return the mean
		virtual double std  ( void ) const { return 1e5; } ///< return the standard deviation
		virtual void shrink ( double xmin, double xmax ) { throw NotImplementedError(); } ///< shrink the prior if it is broader than the range between xmin and xmax
		virtual int get_code(void) const { throw NotImplementedError(); } ///< return the typcode of this prior
		virtual double cdf ( double x ) const { throw NotImplementedError(); } ///< cdf of the prior
		virtual double getprm ( unsigned int prm ) const { throw NotImplementedError(); }
		virtual double ppf ( double p, double start=NULL ) const { throw NotImplementedError(); }
};

/** \brief Uniform prior on an interval
 *
 * This prior defines a uniform distribution on an interval. It's pdf is thus \f$u-l\f$ if \f$x\in(l,u)\f$ and
 * 0 otherwise.
 */
class UniformPrior : public PsiPrior
{
	private:
		double lower;
		double upper;
		double height;
		UniformRandom rng;
	public:
		UniformPrior ( double low, double high ) : lower(low), upper(high), rng ( low, high ) { height = 1./(high-low); } ///< Set up a UniformPrior on the interval from low to high
		UniformPrior ( const UniformPrior& original ) : lower(original.lower),
                                                        upper(original.upper),
                                                        height(original.height),
                                                        rng(original.rng) {} ///< copy constructor
		double pdf ( double x ) const { return ( x>lower && x<upper ? height : 0 ); }                      ///< evaluate the pdf of the prior at position x
		double dpdf ( double x ) { return ( x!=lower && x!=upper ? 0 : (x==lower ? 1e20 : -1e20 ));} ///< derivative of the pdf of the prior at position x (jumps at lower and upper are replaced by large numbers)
		double rand ( void ) { return rng.draw(); }                                                 ///< draw a random number
        PsiPrior * clone ( void ) const { return new UniformPrior(*this); }
		double mean ( void ) const { return 0.5*(lower+upper); }  ///< return the mean
		double std  ( void ) const { return sqrt((upper-lower)*(upper-lower)/12); }
		void shrink ( double xmin, double xmax ) {}         ////< shrinking is not really defined in this case ~> do not shrink
		int get_code(void) const { return 0; } /// return the typcode of this prior
		double cdf ( double x ) const { return ( x<lower ? 0 : (x>upper ? 1 : (x-lower)/(upper-lower) ) ); }
		double getprm ( unsigned int prm ) const { return (prm==0 ? lower : upper ); }
		double ppf ( double p, double start=NULL ) const { return ( p>1 ? upper : (p<0 ? lower : p*(upper-lower)+lower)); }
};

/** \brief gaussian (normal) prior
 *
 * This defines a gaussian prior on the entire real axis. It's pdf is defined by the normal
 *
 \f[
 f(x) = \frac{1}{\sqrt{2 \pi}\sigma} \exp \left( - \frac{(x-\mu)^2}{2\sigma^2} \right).
 \f]
 */
class GaussPrior : public PsiPrior
{
	private:
		double mu;
		double sg;
		double normalization;
		double var;
		double twovar;
		GaussRandom rng;
	public:
		GaussPrior ( double mean, double sd ) : mu(mean), sg(sd), var(sg*sg), twovar(2*sg*sg), rng(mean,sd) { normalization = 1./(sqrt(2*M_PI)*sg); }        ///< initialize prior to have mean mean and standard deviation sd
		GaussPrior ( const GaussPrior& original ) : mu(original.mu),
                                                    sg(original.sg),
                                                    normalization(original.normalization),
                                                    var(original.var),
                                                    twovar(original.twovar),
                                                    rng(original.rng) {} ///< copy contructor
		double pdf ( double x ) const { return normalization * exp ( - (x-mu)*(x-mu)/twovar ); }                                              ///< return pdf of the prior at position x
		double dpdf ( double x ) { return - x * pdf ( x ) / var; }                                                                      ///< return derivative of the prior at position x
		double rand ( void ) {return rng.draw(); }
        PsiPrior * clone ( void ) const { return new GaussPrior(*this); }
		double mean ( void ) const { return mu; } ///< mean
		double std  ( void ) const { return sg; } ///< return standard deviation
		void shrink ( double xmin, double xmax );
		int get_code(void) const { return 1; } /// return the typcode of this prior
		double cdf ( double x ) const { return Phi ( (x-mu)/sg ); }
		double getprm ( unsigned int prm ) const { return ( prm==0 ? mu : sg ); }
		double ppf ( double p, double start=NULL ) const;
};

/** \brief beta prior
 *
 * This defines a beta prior distribution that is defined on an interval. It's pdf is defined by
 *
 \f[
 f(x) = \frac{x^{\alpha-1} (1-x)^{\beta-1}} {B(\alpha,\beta)},
 \f]
 *
 * if \f$x\in(0,1)\f$. It is zero otherwise.
 */
class BetaPrior : public PsiPrior
{
	private:
		double alpha;
		double beta;
		double normalization;
		BetaRandom rng;
		double mode;
	public:
		BetaPrior ( double al, double bt ) : alpha(al), beta(bt), normalization(betaf(al,bt)), rng (alpha, beta) {
            mode = (al-1)/(al+bt-2);
            mode = pdf(mode); }                      ///< Initialize with parameters alpha=al, beta=bt
        BetaPrior ( const BetaPrior& original) : alpha(original.alpha),
                                                 beta(original.beta),
                                                 normalization(original.normalization),
                                                 rng(original.rng),
                                                 mode(original.mode) {} ///< copy constructor
		double pdf ( double x ) const { return (x<0||x>1 ? 0 : pow(x,alpha-1)*pow(1-x,beta-1)/normalization); }             ///< return beta pdf
		double dpdf ( double x ) { return (x<0||x>1 ? 0 : ((alpha-1)*pow(x,alpha-2)*pow(1-x,beta-1) + (beta-1)*pow(1-x,beta-2)*pow(x,alpha-1))/normalization); }      ///< return derivative of beta pdf
		double rand ( void ) {return rng.draw();};                                                                                         ///< draw a random number using rejection sampling
        PsiPrior * clone ( void ) const { return new BetaPrior(*this); }
		double mean ( void ) const { return alpha/(alpha+beta); }
		double std  ( void ) const { return sqrt ( alpha*beta/((alpha+beta)*(alpha+beta)*(alpha+beta+1)) ); }
		void shrink ( double xmin, double xmax );
		int get_code(void) const { return 2; } /// return the typcode of this prior
		double cdf ( double x ) const { return (x<0 ? 0 : (x>1 ? 1 : betainc ( x, alpha, beta ))); }
		double getprm ( unsigned int prm ) const { return ( prm==0 ? alpha : beta ); }
		double ppf ( double p, double start=NULL ) const;
};

/** \brief gamma prior
 *
 * This defines a gamma prior that is defined for the positive axis. It's pdf is defined by
 *
 \f[
 f(x) = x^{k-1} \frac{\exp(-x/\theta)}{\theta^k\Gamma(k)},
 \f]
 * for positive numbers and it is zero otherwise.
 */
class GammaPrior : public PsiPrior
{
	private:
		double k;
		double theta;
		double normalization;
		GammaRandom rng;
	public:
		GammaPrior ( double shape, double scale ) : k(shape), theta(scale), rng(shape, scale) { normalization = exp(gammaln(shape) + shape*log(scale));}                         ///< Initialize a gamma prior
        GammaPrior ( const GammaPrior& original ) : k(original.k),
                                                    theta(original.theta),
                                                    normalization(original.normalization),
                                                    rng(original.rng) {} ///< copy constructor
		virtual double pdf ( double x ) const { return (x>0 ? pow(x,k-1)*exp(-x/theta)/normalization : 0 );}                                                             ///< return pdf at position x
		virtual double dpdf ( double x ) { return (x>0 ? ( (k-1)*pow(x,k-2)*exp(-x/theta)-pow(x,k-1)*exp(-x/theta)/theta)/normalization : 0 ); }                   ///< return derivative of pdf
		virtual double rand ( void ) {return rng.draw(); };
        PsiPrior * clone ( void ) const { return new GammaPrior(*this); }
		virtual double mean ( void ) const { return k*theta; }
		double std  ( void ) const { return sqrt ( k*theta*theta ); }
		void shrink ( double xmin, double xmax );
		virtual int get_code(void) const { return 3; } /// return the typcode of this prior
		virtual double cdf ( double x ) const { return ( x<0 ? 0 : gammainc ( k, x/theta ) / exp ( gammaln ( k ) ) ); }
		double getprm ( unsigned int prm ) const { return ( prm==0 ? k : theta ); }
		virtual double ppf ( double p, double start=NULL ) const;
};

/** \brief negative gamma prior
 *
 * This defines a gamma prior that is defined for the negative axis. Thus, the pdf is defined by
 *
 \f[
 f(x) = (-x)^{k-1} \frac{\exp(x/\theta)}{\theta^k\Gamma(k)},
 \f]
 * for negative numbers and is zero otherwise.
 */
class nGammaPrior : public GammaPrior
{
	public:
		nGammaPrior ( double shape, double scale ) : GammaPrior(shape,scale) {}
        nGammaPrior ( const nGammaPrior& original ) : GammaPrior(original) {} ///< copy constructor
		double pdf ( double x ) const { return GammaPrior::pdf ( -x ); }
		double dpdf ( double x ) { return -GammaPrior::dpdf ( -x ); }
		double rand ( void ) { return -GammaPrior::rand(); }
        PsiPrior * clone ( void ) const { return new nGammaPrior(*this); }
		double mean ( void ) const { return -GammaPrior::mean(); }
		void shrink ( double xmin, double xmax );
		int get_code(void) const { return 4; } /// return the typcode of this prior
		double cdf ( double x ) const { return ( x>0 ? 1 : 1-GammaPrior::cdf ( -x ) ); }
		double ppf ( double p, double start=NULL ) const { return - GammaPrior::ppf ( 1-p ); }
};

/** \brief inverse gamma prior
 *
 * This defines an inverse gamma prior that is conjugate to the variance of a normal distribution
 * and can probably be assumed quasi conjugate to the width or parameters like that. Its pdf is defined
 * by
 *
 \f[
 f(x) = \frac{\beta^\alpha}{\Gamma(\alpha)} x^{-\alpha-1} \exp ( -\beta/x )
 \f]
 */
class invGammaPrior : public PsiPrior
{
	private:
		double alpha;
		double beta;
		double normalization;
		GammaRandom rng;
	public:
		invGammaPrior ( double shape, double scale ) : alpha ( shape ), beta ( scale ), rng ( shape, 1./scale ) { normalization = pow(alpha,beta)/exp(gammaln(shape));} ///< Initialize inverse gamma prior
		invGammaPrior ( const invGammaPrior& original ) : alpha(original.alpha),
					beta(original.beta),
					normalization ( original.normalization ),
					rng(original.rng) {} ///< copy constructor
		virtual double pdf ( double x ) { return ( x>0 ? pow ( x, -alpha-1 ) * exp ( -beta/x ) * normalization : 0 ); }
		virtual double dpdf ( double x ) { return (x>0 ? ( (-alpha-1)*pow(x,-alpha-2) * exp ( -beta/x ) + pow(x,-alpha-1) * exp ( -beta/x ) * beta / (x*x) ) * normalization : 0 ); }
		virtual double rand ( void ) { return 1./rng.draw(); }
		PsiPrior * clone ( void ) const { return new invGammaPrior(*this); }
		virtual double mean ( void ) const { return beta/(alpha-1); }
		double std ( void ) const { return ( alpha>2 ? beta / ( (alpha-1)*sqrt(alpha-2) ) : 1e5 ); }
		virtual void shrink ( double xmin, double xmax ) {} /// Doesn't shrink!!
		virtual int get_code ( void ) const { return 5; } /// return the typecode of this prior
};

/** \brief negative inverse gamma prior
 *
 * This defines a negative inverse gamma prior that is conjugate to the variance of a normal distribution
 * and can probably be assumed quasi conjugate to the width or parameters like that. Its pdf is defined
 * by
 *
 \f[
 f(x) = \frac{\beta^\alpha}{\Gamma(\alpha)} (-x)^{-\alpha-1} \exp ( \beta/x )
 \f]
 */
class ninvGammaPrior : public invGammaPrior
{
	public:
		ninvGammaPrior ( double shape, double scale ) : invGammaPrior ( shape, scale ) {}
		ninvGammaPrior ( const ninvGammaPrior& original ) : invGammaPrior ( original ) {}
		double pdf ( double x ) { return invGammaPrior::pdf ( -x ); }
		double dpdf ( double x ) { return -invGammaPrior::dpdf ( -x ); }
		double rand ( void ) { return -invGammaPrior::rand(); }
		PsiPrior * clone ( void ) const { return new ninvGammaPrior ( *this ); }
		double mean ( void ) const { return -invGammaPrior::mean(); }
		void shrink ( double xmin, double xmax ) { invGammaPrior::shrink ( -xmax, -xmin ); }
		int get_code ( void ) const { return 6; } /// return the typecode of this prior
};
#endif
