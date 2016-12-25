/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#ifndef SIGMOID_H
#define SIGMOID_H

#include "errors.h"
#include "special.h"
#include <cmath>
#include <string>

/** \brief common base class for all sigmoids */
class PsiSigmoid
{
	public:
		virtual double f   ( double x ) const { throw NotImplementedError(); }            ///< This should return the value of the sigmoid itself (between 0 and 1)
		virtual double df  ( double x ) const { throw NotImplementedError(); }            ///< This should give the first derivative of the sigmoid
		virtual double ddf ( double x ) const { throw NotImplementedError(); }            ///< This should give the second derivative of the sigmoid
		virtual double inv ( double p ) const { throw NotImplementedError(); }            ///< This should give the inverse of the sigmoid (taking values between 0 and 1)
		virtual int    getcode ( void ) const { throw NotImplementedError(); }            ///< return the sigmoid identifier
        virtual PsiSigmoid * clone ( void ) const { throw NotImplementedError(); }                ///< clone object by value
        static std::string getDescriptor ( void ) { throw NotImplementedError(); }///< get a short string that identifies the type of sigmoid
};

/** \brief identity function as sigmoid
 *
 * This is useful if the desired function cannot be separated to the form sigmoid + core
 * and the whole function has to be implemented in the core object. An example for this
 * case is the Naka-Rushton nonlinearity typically used to model electrophysiological
 * data.
 */
class PsiId : public PsiSigmoid
{
	public:
		double f   ( double x ) const { return x; }
		double df  ( double x ) const { return 1; }
		double ddf ( double x ) const { return 0; }
		double inv ( double x ) const { return x; }
		int getcode ( void ) const { return 6; }
		PsiSigmoid * clone ( void ) const { return new PsiId(*this); }
		static std::string getDescriptor ( void ) { return "id"; }
};

/** \brief logistic function
 *
 * The logistic function is given by f(x) = 1/(1+exp(-x))
 */
class PsiLogistic : public PsiSigmoid
{
	public:
		PsiLogistic ( void ) {}  ///< constructor
		PsiLogistic ( const PsiLogistic& original) {}  ///< copy constructor
		double f ( double x ) const;                 ///< value of the sigmoid at position x
		double df ( double x ) const;                ///< derivative of the sigmoid at position x
		double ddf ( double x ) const;               ///< second derivative of the sigmoid
		double inv ( double p ) const { return log(p/(1-p)); }  ///< inverse of the sigmoid
		int getcode ( void ) const { return 1; }     ///< return the sigmoid identifier
        PsiSigmoid * clone ( void ) const {
            return new PsiLogistic(*this);
        }
        static std::string getDescriptor ( void ) {
            return "logistic";
        }
};

/** \brief gaussian cdf function
 *
 * The gaussian cdf function is given by f(x) = Phi(x), where Phi is the cumulative distribution function for the gaussian distribution.
 */
class PsiGauss : public PsiSigmoid
{
	public:
		PsiGauss ( void ) {} ///< constructor
		PsiGauss ( const PsiGauss& original) {} ///< copy constructor
		double f   ( double x ) const;                 ///< value of the sigmoid at x
		double df  ( double x ) const;                 ///< derivative of the sigmoid at x
		double ddf ( double x ) const;                 ///< second derivative of the sigmoid at x
		double inv ( double p ) const;                 ///< inverse of the sigmoid
		int getcode ( void ) const { return 2; }       ///< return the sigmoid identifier
        PsiSigmoid * clone (void ) const {
            return new PsiGauss(*this);
        }
        static std::string getDescriptor ( void ) {
            return "gauss";
        }
};

/** \brief left-skewed gumbel cdf
 *
 * Cumulative densitiy function of the gumbel distribution.
 */
class PsiGumbelL : public PsiSigmoid
{
	public:
		PsiGumbelL ( void ) {} ///< contructor
		PsiGumbelL ( const PsiGumbelL& original ) {} ///< copy constructor
		double f   ( double x ) const;              ///< returns the value of the gumbel cdf at position x
		double df  ( double x ) const;              ///< returns the derivative of the gumbel cdf at position x
		double ddf ( double x ) const;              ///< returns the 2nd derivative of the gumbel cdf at position x
		double inv ( double p ) const;              ///< returns the inverse of the gumbel cdf at position p
		int getcode ( void ) const { return 3; }    ///< return the sigmoid identifier
        PsiSigmoid * clone ( void ) const {
            return new PsiGumbelL(*this);
        }
        static std::string getDescriptor ( void ) {
            return "gumbel_l";
        }
};

/** \brief right-skewed gumbel cdf
 *
 * Cumulative densitiy function of the gumbel distribution.
 */
class PsiGumbelR : public PsiSigmoid
{
	public:
		PsiGumbelR ( void ) {} ///< constructor
		PsiGumbelR ( const PsiGumbelR& original ) {} ///< copy constructor
		double f   ( double x ) const;             ///< returns the value of the right skewed gumbel cdf at position x
		double df  ( double x ) const;             ///< returns the derivative of the right skewed gumbel cdf at position x
		double ddf ( double x ) const;             ///< returns the 2nd derivative of the right skewed gumbel cdf at position x
		double inv ( double p ) const;             ///< returns the inverse of the right skewed gumbel cdf at position p
		int getcode ( void ) const { return 6; }   ///< return the sigmoid identifier
        PsiSigmoid * clone ( void ) const {
            return new PsiGumbelR(*this);
        }
        static std::string getDescriptor ( void ) {
            return "gumbel_r";
        }
};

/** \brief cauchy cdf
 *
 * Cumulative density function of the cauchy distribution
 */
class PsiCauchy : public PsiSigmoid
{
	public:
        PsiCauchy( void ) {}                 ///< constructor
        PsiCauchy( const PsiCauchy& oiginal) {} ///< copy constructor
		double f   ( double x ) const;             ///< returns the value of the cauchy cdf at position x
		double df  ( double x ) const;             ///< returns the derivative of the cauchy cdf at position x
		double ddf ( double x ) const;             ///< returns the 2nd derivative of the cauchy cdf at position x
		double inv ( double p ) const;             ///< returns the inverse of the cauchy cdf at position x
		int    getcode ( void ) const { return 4; }///< returns the sigmoid identifier
        PsiSigmoid * clone ( void ) const {
            return new PsiCauchy(*this);
        }
        static std::string getDescriptor ( void ) {
            return "cauchy";
        }
};

/** \brief exponential cdf
 *
 * Cumulative density function of the exponential distribution
 * combined with a polyCore this will give a weibull
 */
class PsiExponential : public PsiSigmoid
{
	public:
        PsiExponential( void ) {}                 ///< constructor
        PsiExponential( const PsiExponential& oiginal) {} ///< copy constructor
		double f   (double x ) const;              ///< returns the value of the exponential cdf at position x
		double df  (double x ) const;              ///< returns the derivative of the exponential cdf at position x
		double ddf (double x ) const;              ///< returns the 2nd derivative of the exponential cdf at position x
		double inv (double p ) const throw(BadArgumentError);              ///< returns the return the inverse of the exponential cdf at position x
		int    getcode ( void ) const { return 5; }///< returns the sigmoid identifier
        PsiSigmoid * clone ( void ) const {
            return new PsiExponential(*this);
        }
        static std::string getDescriptor ( void ) {
            return "exponential";
        }
};

#endif
