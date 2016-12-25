/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#ifndef PSYCHOMETRIC_H
#define PSYCHOMETRIC_H

#include <vector>
#include <cmath>
#include <algorithm>
#include "core.h"
#include "sigmoid.h"
#include "errors.h"
#include "prior.h"
#include "data.h"
#include "linalg.h"

/** \brief Standard psychometric function model
 *
 * Standard model for the psychometric function that assumes that the number of correct responses is a
 * binomial random variable with parameters N (number of trials) and Psi, where Psi is
 *
 * Psi = guessingrate + (1-guessingrate-lapserate) * Sigmoid ( x | theta )
 *
 * For an nAFC task, the guessingrate is typicall fixed at 1/n.
 *
 * The term Sigmoid ( x | theta ) is represented by two objects:
 * a PsiSigmoid f, that describes a nonlinear function from the real numbers to (0,1) and
 * a PsiCore g, that describes the "internal" workings of the nonlinear function. Thus,
 * the Term Sigmoid ( x | theta ) equals f(g(x,theta)).
 *
 * The object copies and manages its own memory. This means that when setting
 * PsiSigmoid and PsiCore in the constructor, as well as setting PsiPrior via
 * setPrior(), copies of the respective objects are made. Furthermore when the
 * object is destroyed theses copies are destoryed also.
 */
class PsiPsychometric {
	private:
		int Nalternatives;
		double guessingrate;
		bool gammaislambda;
		PsiCore * Core;
		PsiSigmoid * Sigmoid;
		std::vector<PsiPrior*> priors;
	protected:
		PsiPsychometric (
			int nAFC,                                                               ///< number of alternatives (1 indicating yes/no)
			PsiCore * core,                                                         ///< internal part of the nonlinear function
			PsiSigmoid * sigmoid,                                                   ///< "external" saturating part of the nonlinear function
			unsigned int nparameters                                                ///< number of parameters given explicitely
			);                  ///< Set up a psychometric function model for an nAFC task, explicitely specifiing the number of parameters (useful for derived classes)
	public:
		PsiPsychometric (
			int nAFC,                                                                ///< number of alternatives in the task (1 indicating yes/no)
			PsiCore * core,                                                          ///< internal part of the nonlinear function (in many cases this is actually a linear function)
			PsiSigmoid * sigmoid                                                     ///< "external" saturating part of the nonlinear function
			);    ///< Set up a psychometric function model for an nAFC task (nAFC=1 ~> yes/no)
		~PsiPsychometric ( void );   ///< destructor (also deletes the core and sigmoid objects)
		virtual double evaluate (
			double x,                                                                ///< stimulus intensity
			const std::vector<double>& prm                                           ///< parameters of the psychometric function model
			) const;  ///< Evaluate the psychometric function at this position
		virtual double negllikeli (
			const std::vector<double>& prm,                                          ///< parameters of the psychometric function model
			const PsiData* data                                                      ///< data for which the likelihood should be evaluated
			) const;   ///< negative log likelihood
		virtual double neglpost (
			const std::vector<double>& prm,                                          ///< parameters of the psychometric function model
			const PsiData* data                                                      ///< data for which the posterior should be evaluated
			) const;     ///< negative log posterior  (unnormalized)
		virtual double leastfavourable (
			const std::vector<double>& prm,                                          ///< parameters of the psychometric function model
			const PsiData* data,                                                     ///< data for which the likelihood should be evaluated
			double cut,                                                              ///< performance level at which the threshold should be evaluated
			bool threshold=true                                                      ///< should the calculations be performed for thresholds? (anything else is not yet implemented)
			) const; ///< derivative of log likelihood in the least favourable direction in parameter space
		virtual double deviance (
			const std::vector<double>& prm,                                          ///< parameters of the psychometric functin model
			const PsiData* data                                                      ///< data for which the likelihood should be evaluated
			) const; ///< deviance for a given data set and parameter constellation
		virtual Matrix * ddnegllikeli (
				const std::vector<double>& prm,                                      ///< parameters at which the second derivative should be evaluated
				const PsiData* data                                                  ///< data for which the likelihood should be evaluated
				) const;                                          ///< 2nd derivative of the negative log likelihood (newly allocated matrix)
		virtual std::vector<double> dnegllikeli (
				const std::vector<double>& prm,                                      ///< parameters at which the first derivative should be evaluated
				const PsiData* data                                                  ///< data for which the likelihood should be evaluated
				) const;                                          ///< 1st derivative of the negative log likelihood
		const PsiCore* getCore ( void ) const { return Core; }                ///< get the core of the psychometric function
		const PsiSigmoid* getSigmoid ( void ) const { return Sigmoid; }       ///< get the sigmoid of the psychometric function
		void setPrior ( unsigned int index, PsiPrior* prior ) throw(BadArgumentError);                   ///< set a Prior for the parameter indicated by index
		double evalPrior ( unsigned int index, double x ) const {return priors[index]->pdf(x);}              ///< evaluate the respective prior at value x
		virtual double randPrior ( unsigned int index ) const { return priors[index]->rand(); }                            ///< sample form a prior
		const PsiPrior* getPrior ( unsigned int index ) const { return priors[index]; } ///< get a prior
		int getNalternatives ( void ) const { return Nalternatives; }         ///< get the number of alternatives (1 means yes/no)
		virtual unsigned int getNparams ( void ) const { return (Nalternatives==1 ? (gammaislambda ? 3 : 4 ) : 3 ); } ///< get the number of free parameters of the psychometric function
		virtual std::vector<double> getStart ( const PsiData* data ) const ;                ///< determine a starting value using logistic regression on a dataset
		double getThres (
			const std::vector<double>& prm,                                          ///< parameters of the psychometric function model
			double cut                                                               ///< performance level at which the threshold should be evaluated
			) const { return Core->inv(Sigmoid->inv(cut),prm); }  ///< get the threshold at a cut between 0 and 1
		double getSlope (
			const std::vector<double>& prm,                                          ///< parameters of the psychometric function model
			double x                                                                 ///< performance level at which the slope should be evaluated
			) const { return Sigmoid->df ( Core->g ( x, prm ) ) * Core->dgx ( x, prm ); } ///< get the slope at a stimulus intensity
		std::vector<double> getDevianceResiduals (
			const std::vector<double>& prm,                                          ///< parameters of the psychometric function model
			const PsiData* data                                                      ///< data for which the deviance residuals should be determined
			) const;  ///< deviance residuals for model checking
		double getRpd (
			const std::vector<double>& devianceresiduals,                            ///< deviance residuals as determined by getDevianceResiduals
			const std::vector<double>& prm,                                          ///< parameters of the psychometric function model
			const PsiData* data                                                      ///< data set corresponding to the deviance residuals
			) const;          ///< correlation between deviance residuals and predictions
		double getRkd ( const std::vector<double>& devianceresiduals, const PsiData* data ) const;        ///< correlation between deviance residuals and block sequence
		double dllikeli (
			std::vector<double> prm,                                                     ///< parameters of the model
			const PsiData* data,                                                         ///< data for which the likelihood should be evaluated
			unsigned int i                                                               ///< index of the parameter for which the derivative should be evaluated
			) const;                                                                 ///< derivative of the negative loglikelihood with respect to parameter i
		double dlposteri (
			std::vector<double> prm,                                                     ///< parameters of the psychometric function model
			const PsiData* data,                                                         ///< data for which the likelihood should be valuated
			unsigned int i                                                               ///< index of the parameter for which the derivative should be evaluated
			) const;                                                                 ///< derivative of the negative log posterior with respect to parameter i
		void setgammatolambda ( void ) { gammaislambda=true; };                          ///< calling this function applies the constraint that gamma and lambda should be equal in a yes/no paradigm
		double getGuess ( const std::vector<double>& prm ) const { return (gammaislambda ? prm[2] : ( getNalternatives() < 2 ? prm[3] : 1./Nalternatives )); }
		double dpredict ( const std::vector<double>& prm, double x, unsigned int i ) const;    ///< partial derivative of psychometric function prediction w.r.t. i-th parameter
		double ddpredict ( const std::vector<double>& prm, double x, unsigned int i, unsigned int j ) const;    ///< 2nd partial derivative of psychometric function prediction w.r.t. i-th and j-th parameters
};

/** \brief Psychometric function that allows for models the variance of the data by a beta distribution
 *
 * When fitting psychometric functions, we usually assume independence of the responses. Is this independence is violated,
 * the binomial variance model is no longer adequate. The beta psychometric function can better deal with this. It fits an
 * additional parameter nu, such that the subjects responses are distributed according to a beta distribution with parameters
 * alpha = Psi(x)*nu*n+1, beta = (Psi(x)-1)*nu*n+1
 */
class BetaPsychometric : public PsiPsychometric {
	private:
		double fznull ( unsigned int z, const PsiData* data, double nu ) const;
		double negllikelinull ( const PsiData* data, double nu ) const;
	public:
		BetaPsychometric ( int nAFC, PsiCore * core, PsiSigmoid * sigmoid ) : PsiPsychometric ( nAFC, core, sigmoid, ( nAFC<2 ? 5 : 4 ) ) {}
		double negllikeli (
			const std::vector<double>& prm,           ///< parameters of the psychometric function model
			const PsiData* data                       ///< data for which the likelihood should be evaluated
			) const; ///< negative log likelihood
		std::vector<double> dnegllikeli (
				const std::vector<double>& prm,       ///< parameters at which the first derivative should be evaluated
				const PsiData* data                   ///< data for which the likelihood should be evaluated
				) const;                 ///< 1st derivative of the negative log likelihood
		Matrix * ddnegllikeli (
				const std::vector<double>& prm,       ///< parameters at which the second derivative should be evaluated
				const PsiData* data                   ///< data for which the likelihood should be evaluated
				) const;                 ///< 2nd derivative of the negative log likelihood (newly allocated matrix)
		unsigned int getNparams ( void ) const { return PsiPsychometric::getNparams()+1; }   ///< get the number of free parameters of the psychometric function
		double deviance (
			const std::vector<double>& prm,                      ///< parameters of the psychometric function model
			const PsiData * data                                 ///< data for which the likelihood should be evaluated
			) const; ///< deviance for a given data set and parameter constellation

		std::vector<double> getStart ( const PsiData* data ) const { std::vector<double> out (PsiPsychometric::getStart ( data )); out[out.size()-1] = .99999; return out;}
};

/** \brief Psychometric function with one separate data point
 *
 * Wichmann & Hill (2001) suggest to detect whether or not a data point x0 is an outlier by
 * fitting a new model to all data points except for the data point of interest and to use
 * a separate parameter to fit x0. If p is the number of correct trials at position x0, that is
 *
 * Psi(x0) = p
 * Psi(x)  = guessingrate + (1-guessingrate-lapserate) * Sigmoid ( x | theta ), if x != x0
 *
 * This model is then compared to the standard psychometric function model.
 *
 * Due to the discontinuity at x0, derivatives in this model are treated with respect to the
 * fitted function.
 */

class OutlierModel : public PsiPsychometric {
	private:
		unsigned int jout;
		double getp ( const std::vector<double>& prm ) const;
	public:
		OutlierModel (
			int nAFC,                                                            ///< number of alternatives in the task (1 indicating yes/no)
			PsiCore * core,                                                      ///< internal part of the nonlinear function
			PsiSigmoid * sigmoid,                                                ///< "external" saturating part of the nonlinear function
			unsigned int exclude                                                 ///< index of the data block to be excluded
			) : PsiPsychometric ( nAFC, core, sigmoid ), jout(exclude) {}; ///< set up a psychometric function model that treats one block separately
		void setexclude ( unsigned int exclude ) { jout = exclude; }   ///< change the excluded block
		double negllikeli (
			const std::vector<double>& prm,                                      ///< parameters of the psychometric function model
			const PsiData * data                                                 ///< data for which the likelihood should be evaluated
			) const;                         ///< negative log likelihood
		double neglpost (
			const std::vector<double>& prm,                                      ///< parameters of the psychometric function model
			const PsiData * data                                                 ///< data for which the likelihood should be evaluated
			) const;                         ///< negative log likelihood
		double deviance (
			const std::vector<double>& prm,                                      ///< parameters of the psychometric function model
			const PsiData* data                                                  ///< data for which the deviance should be evaluated
			) const;                        ///< deviance
		unsigned int getNparams ( void ) const { return PsiPsychometric::getNparams()+1; }
		double randPrior ( unsigned int index ) const { return ( index<PsiPsychometric::getNparams() ? PsiPsychometric::randPrior(index) : PsiRandom().rngcall() ); }                            ///< sample form a prior
};

#endif
