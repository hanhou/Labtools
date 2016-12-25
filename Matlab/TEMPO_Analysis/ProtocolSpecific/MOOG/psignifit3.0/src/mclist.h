/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#ifndef MCLIST_H
#define MCLIST_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include "errors.h"
#include "special.h"
#include "data.h"
#include "rng.h"

/** \brief basic monte carlo samples list
 *
 * This list stores monte carlo samples and deviances, nothing else.
 */
class PsiMClist
{
	private:
		std::vector< std::vector<double> > mcestimates;
		std::vector<double> deviances;
	public:
		PsiMClist (
			int N,                      ///< number of samples to be drawn
			int nprm                    ///< number of parameters in the model that is analyzed
			) : mcestimates(nprm, std::vector<double>(N) ), deviances(N) {}   ///< Initialize the list to take N samples of nprm parameters
		PsiMClist ( const PsiMClist& mclist ) : mcestimates ( mclist.mcestimates ), deviances ( mclist.deviances ) {}   ///< copy a list of mcsamples
		~PsiMClist ( ) {} ///< destructor
		std::vector<double> getEst ( unsigned int i ) const;       ///< get a single parameter estimate at sample i
		double getEst (
			unsigned int i,                                        ///< sample index
			unsigned int prm                                       ///< parameter index
			) const;                                                           ///< get a single sample of a single parameter
		void setEst (
			unsigned int i,                                        ///< index of the sample to be set
			const std::vector<double> est,                ///< parameter vector to be set at index
			double deviance                               ///< deviance associated with the sample
			);                                                                 ///< set a sample of parameters
		virtual void setdeviance ( unsigned int i, double deviance );                   ///< set the deviance separately for sample i
		virtual double getPercentile (
			double p,                                     ///< desired percentile (in the range (0,1))
			unsigned int prm                                       ///< index of the paramter of interest
			);                                                                 ///< get a percentile for parameter prm
		virtual double getMean (
			unsigned int prm                              ///< index of the parameter of interest
			) const ;                                                          ///< get the average of parameter prm
		virtual double getStd (
			unsigned int prm                             ///< index of the parameter of interest
			) const ;                                                          ///< get the standard deviantion of parameter prm
		double getdeviance ( unsigned int i ) const;                                    ///< get the deviance of sample i
		unsigned int getNsamples ( void ) const { return mcestimates[0].size(); }       ///< get the total number of samples
		unsigned int getNparams ( void ) const { return mcestimates.size(); }           ///< get the number of parameters
		double getDeviancePercentile ( double p );                             ///< get the p-percentile of the deviance (p in the range (0,1) )
};

/** \brief list of bootstrap samples
 *
 * Bootstrap samples support some special operations that regular monte carlo samples don't. In particular bootstrap
 * samples from a psychometric function should incorporate information about
 *
 * -# The thresholds that are associated with each parameter vector
 * -# the bootstrap samples themselves and not only the resulting parameter estimates
 * -# correlations of the psychometric function with the bootstrap samples "sequence"
 */
class BootstrapList : public PsiMClist
{
	private:
		bool BCa;
		std::vector<double> acceleration_t;
		std::vector<double> bias_t;
		std::vector<double> acceleration_s;
		std::vector<double> bias_s;
		std::vector< std::vector<int> > data;
		std::vector<double> cuts;
		std::vector< std::vector<double> > thresholds;
		std::vector< std::vector<double> > slopes;
		std::vector<double> Rpd;
		std::vector<double> Rkd;
	public:
		BootstrapList (
			unsigned int N,                                              ///< number of samples to be drawn
			unsigned int nprm,                                           ///< number of parameters in the model
			unsigned int nblocks,                                        ///< number of blocks in the experiment
			std::vector<double> Cuts                                     ///< performance levels at which thresholds should be determined
			) : PsiMClist (N,nprm),
				BCa(false),
				acceleration_t(Cuts.size()),
				bias_t(Cuts.size()),
				acceleration_s(Cuts.size()),
				bias_s(Cuts.size()),
				data(N,std::vector<int>(nblocks)),
				cuts(Cuts),
				thresholds (Cuts.size(), std::vector<double> (N)),
				slopes     (Cuts.size(), std::vector<double> (N)),
				Rpd(N),
				Rkd(N)
			{ } ///< set up the list
		// TODO: should setBCa be private and friend of parametric bootstrap?
		void setBCa_t (
			unsigned int i,                                               ///< index of the cut for which Bias and Acceleration should be set
			double Bias,                                                  ///< Bias to be set
			double Acceleration                                           ///< Acceleration to be set
			) { BCa=true; bias_t[i] = Bias; acceleration_t[i] = Acceleration; }  ///< set bias and acceleration to get BCa confidence intervals
		void setBCa_s (
			unsigned int i,                                               ///< index of the cut for which Bias and Acceleration should be set
			double Bias,                                                  ///< Bias to be set
			double Acceleration                                           ///< Acceleration to be set
			) { BCa=true; bias_s[i] = Bias; acceleration_s[i] = Acceleration; }  ///< set bias and acceleration to get BCa confidence intervals
		void setData (
			unsigned int i,                                               ///< index of the bootstrap sample to be set
			const std::vector<int>& newdata                                ///< response counts in the new bootstrap sample (not proportion correct)
			);   ///< store a simulated data set
		std::vector<int> getData ( unsigned int i ) const;                 ///< get a simulated data set at posititon i

		double getThres ( double p, unsigned int cut );                    ///< get the p-th percentile associated with the threshold at cut
		double getThres_byPos ( unsigned int i, unsigned int cut );        ///< get the threshold for the i-th sample
		void setThres ( double thres, ///< new value of the threshold
                unsigned int i,       ///< index of the bootstrap sample
                unsigned int cut      ///< index of the desired cut
                );  ///< set the value of a threshold associated with the threshold at cut

		double getSlope ( double p, unsigned int cut );                    ///< get the p-th percentile associated with the slope at the given cut
		double getSlope_byPos ( unsigned int i, unsigned int cut );        ///< get the slope at cut for the i-th sample
		void setSlope ( double sl,    ///< new value of the slope
				unsigned int i,       ///< index of the bootstrap sample
				unsigned int cut      ///< index of the desired cut
				);  ///< set the value of the slope associated with the threshold at cut

		unsigned int getNblocks ( void ) const { return data[0].size(); }  ///< get the number of blocks in the underlying dataset
		double getCut ( unsigned int i ) const;                            ///< get the value of cut i
		double getAcc_t ( unsigned int i ) const { return acceleration_t[i]; } ///< get the acceleration constant for cut i
		double getBias_t ( unsigned int i ) const { return bias_t[i]; }       ///< get the bias for cut i
		double getAcc_s ( unsigned int i ) const { return acceleration_s[i]; } ///< get the acceleration constant for cut i
		double getBias_s ( unsigned int i ) const { return bias_s[i]; }       ///< get the bias for cut i
		// TODO: should setRpd be private and friend of parametricbootstrap?
		void setRpd ( unsigned int i, ///< index of the bootstrap sample
                double r_pd           ///< new correlation
                );///< set correlation between predicted values and deviance residuals for a simulated dataset
		double getRpd ( unsigned int i ) const;                            ///< get correlation between predicted values and deviance residuals for simulated dataset i
		double percRpd ( double p );                                       ///< get the p-th percentile of the correlations between predicted values and deviance residuals
		// TODO: should setRkd be private and friend of parametric bootstrap?
		void setRkd ( unsigned int i, double r_kd );                       ///< set correlation between block index and deviance residuals for a simulated dataset
		double getRkd ( unsigned int i ) const;                            ///< get correlation between block index and deviance residuals for simulated dataset i
		double percRkd ( double p );                                       ///< get the p-th percentile of the correlations between block index and deviance residuals
};

/** \brief list of JackKnife data
 *
 * JackKnifeing is not suggested for the assessment of confidence intervals or variability. Instead the close link between jackknife samples
 * and individual data points is useful to determine influential data points and outliers.
 */
class JackKnifeList : public PsiMClist
{
	private:
		double maxdeviance;
		std::vector<double> mlestimate;
	public:
		JackKnifeList (
			unsigned int nblocks,                                             ///< number of blocks in the experiment
			unsigned int nprm,                                                ///< number of parameters in the model
			double maxldev,                                                   ///< deviance of the maximum likelihood estimate on the full dataset
			std::vector<double> maxlest                                       ///< maximum likelihood estimate of the full dataset
			) : PsiMClist ( nblocks, nprm ), maxdeviance(maxldev), mlestimate(maxlest) {}    ///< constructor
		unsigned int getNblocks ( void ) const { return getNsamples(); } ///< get the number of blocks in the current experiment
		/** determination of influential observations is performed by checking whether a parameter changes significantly (as defined by
		 * the confidence intervals) if one observation is omitted. Thus, if leaving out one observation results in significant changes
		 * in the estimated parameters, this observation is considered "influential".
		 *
		 * \param block     index of the block to be checked
		 * \param estimate  point estimate of the parameters in the model
		 * \param ci_lower  lower confidence limits for each parameter in the model
		 * \param ci_upper  upper confidence limits for each parameter in the model
		 *
		 * \return a number indicating the influence of the block. Values > 1 correspond to point estimates for that block that are precisely on the CI limits
		 */
		double influential ( unsigned int block, const std::vector<double>& ci_lower, const std::vector<double>& ci_upper ) const;
		/** determination of outliers is based on the following idea: We add a new parameter that fits the data in block perfectly.
		 * If this "modified" model is significantly better than the original model, then this block is considered an outlier.
		 *
		 * \param block      index of the block to be checked
		 *
		 * \return true if block presents an outlier
		 */
		bool outlier ( unsigned int block ) const ; ///< is block an outlier?
};

/** \brief a list of Bayesian MCMC samples
 *
 * This list stores additional data that are important for bayesian analysis.:
 * 1. For each parameter sample, a sample from the respective psychometric function is stored. These samples
 *    can be considered samples from the posterior predictive distribution
 * 2. For each sample from the posterior predictive distribution, the deviance is stored
 * 3. the list allows to obtain the estimated bayesian p-value
 */
class MCMCList : public PsiMClist
{
	private:
		std::vector<double> posterior_Rpd;
		std::vector<double> posterior_Rkd;
		std::vector< std::vector<int> > posterior_predictive_data;
		std::vector<double> posterior_predictive_deviances;
		std::vector<double> posterior_predictive_Rpd;
		std::vector<double> posterior_predictive_Rkd;
		std::vector< std::vector<double> > logratios;       // log ratios of the unnormalized posteriors for the full model and the models with one block omitted
		double accept_rate;
		double H;
	public:
		MCMCList (
			unsigned int N,                                                ///< number of samples to be drawn
			unsigned int nprm,                                             ///< number of parameters in the model
			unsigned int nblocks                                           ///< number of blocks in the experiment
			) : PsiMClist ( N, nprm),
				posterior_Rpd(N),
				posterior_Rkd(N),
				posterior_predictive_data(N,std::vector<int>(nblocks)),
				posterior_predictive_deviances ( N ),
				posterior_predictive_Rpd ( N ),
				posterior_predictive_Rkd ( N ),
				logratios ( N, std::vector<double>(nblocks ) ) {};      ///< set up MCMCList
		void setppData (
			unsigned int i,                                                ///< index of the posterior predictive sample to be set
			const std::vector<int>& ppdata,                                ///< posterior predictive data sample
			double ppdeviance                                              ///< deviance associated with the posterior predictive sample
			);               ///< store a posterior predictive data set
		std::vector<int> getppData ( unsigned int i ) const;               ///< get a posterior predictive data sample
		int getppData ( unsigned int i, unsigned int j ) const;
		double getppDeviance ( unsigned int i ) const;                     ///< get deviance associated with a posterior predictive sample
		void setppRpd ( unsigned int i, double Rpd );
		double getppRpd ( unsigned int i ) const;
		void setppRkd ( unsigned int i, double Rkd );
		double getppRkd ( unsigned int i ) const;
		void setRpd ( unsigned int i, double Rpd );
		double getRpd ( unsigned int i ) const;
		void setRkd ( unsigned int i, double Rkd );
		double getRkd ( unsigned int i ) const;
		unsigned int getNblocks ( void ) const { return posterior_predictive_data[0].size(); }  ///< get the number of blocks
		void setlogratio ( unsigned int i, unsigned int j, double logratio );              ///< set the log posterior ratio for sample i and block j
		double getlogratio ( unsigned int i, unsigned int j ) const;                       ///< get the log posterior ratio for sample i and block j
		void set_accept_rate(double rate) {accept_rate = rate; }  ///< set the acceptance rate
		double get_accept_rate(void) const {return accept_rate; } ///< get the acceptance rate
		void set_entropy ( double entropy ) { H = entropy; } ///< set the entropy if needed
		double get_entropy ( void ) const { return H; }
};

void newsample ( const PsiData * data, const std::vector<double>& p, std::vector<int> * sample );

#endif
