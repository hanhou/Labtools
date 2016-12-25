/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#ifndef MCMC_H
#define MCMC_H

#include <vector>
#include "psychometric.h"
#include "rng.h"
#include "mclist.h"
#include "getstart.h"

class PsiSampler
{
	private:
		const PsiPsychometric * model;
		const PsiData * data;
	public:
		PsiSampler ( const PsiPsychometric * Model, const PsiData * Data ) : model(Model), data(Data) {}///< set up a sampler to sample from the posterior of the parameters of pmf given the data dat
		virtual std::vector<double> draw ( void ) { throw NotImplementedError(); }                     ///< draw a sample from the posterior
		virtual void setTheta ( const std::vector<double> theta ) { throw NotImplementedError(); }     ///< set the "state" of the underlying markov chain
		virtual std::vector<double> getTheta ( void ) { throw NotImplementedError(); }                 ///< get the "state" of the underlying markov chain
		virtual void setStepSize ( double size, unsigned int param ) { throw NotImplementedError(); }  ///< set the size of the steps for parameter param of the sampler
		virtual void setStepSize ( const std::vector<double>& sizes ) { throw NotImplementedError(); } ///< set all stepsizes of the sampler
		virtual double getDeviance ( void ) { throw NotImplementedError(); }                           ///< return the model deviance for the current state
		virtual MCMCList sample ( unsigned int N ) { throw NotImplementedError(); }                   ///< draw N samples from the posterior
		const PsiPsychometric * getModel() const { return model; }                                     ///< return the underlying model instance
		const PsiData         * getData()  const { return data;  }                                     ///< return the underlying data instance
};

class MetropolisHastings : public PsiSampler
{
	private:
		PsiRandom* propose;
		std::vector<double> currenttheta;
		std::vector<double> newtheta;
		std::vector<double> stepwidths;
		double currentdeviance;
		int accept;
	protected:
		double qold;
	public:
		MetropolisHastings (
			const PsiPsychometric * Model,                                                  ///< psychometric funciton model to sample from
			const PsiData * Data,                                                           ///< data to base inference on
			PsiRandom* proposal                                                             ///< proposal distribution (will usually be a gaussian)
			);                                                          ///< initialize the sampler
		~MetropolisHastings ( void ) { delete propose; }
		std::vector<double> draw ( void );                                                ///< perform a metropolis hastings step and draw a sample from the posterior
		virtual double acceptance_probability ( const std::vector<double>& current_theta, const std::vector<double>& new_theta );
		void setTheta ( const std::vector<double>& prm );                                 ///< set the current state of the sampler
		std::vector<double> getTheta ( void ) { return currenttheta; }                    ///< get the current state of the sampler
		double getDeviance ( void ) { return currentdeviance; }                           ///< get the current deviance
		void setStepSize ( double size, unsigned int param );                             ///< set the standard deviation of the proposal distribution for parameter param
		void setStepSize ( const std::vector<double>& sizes );                            ///< set standard deviations of the proposal distribution for all parameters at once
		std::vector<double> getStepsize ( void ) { return stepwidths; }		  			  ///< return the current stepwidth (standard deviations of the proposal distribution)
		MCMCList sample ( unsigned int N );                                               ///< draw N samples from the posterior
		unsigned int getNparams ( void ) { return newtheta.size(); }                      ///< get the number of parameters for which the sampler is set up
		virtual void proposePoint( std::vector<double> &current_theta,
									std::vector<double> &step_widths,
									PsiRandom * proposal,
									std::vector<double> &new_theta);		  			  ///< propose a new sample and save it in new_theta
};

class GenericMetropolis : public MetropolisHastings
{
	private:
		int currentindex;													  			  ///< parameter index for proposing a new sample
	public:
		GenericMetropolis (
			const PsiPsychometric * Model,                                    			  ///< psychometric funciton model to sample from
			const PsiData * Data,                                             			  ///< data to base inference on
			PsiRandom* proposal                                               			  ///< proposal distribution (will usually be a gaussian)
			): MetropolisHastings ( Model, Data, proposal ),
			   currentindex(0) {}
		void proposePoint( std::vector<double> &current_theta,
							std::vector<double> &step_widths,
							PsiRandom * proposal,
							std::vector<double> &new_theta);				  			  ///< propose a new sample and save it in new_theta
		/** \brief Find the optimal stepwidth by regressing each parameter against the others.
		 *
		 * For each parameter, do a regression using QR-decomposition and
		 * take the residuals to calculate the optimal stepwidth.
		 *
		 * @param pilot a pilot sample to base the regression on
		 */
		void findOptimalStepwidth ( PsiMClist const &pilot );
};

class DefaultMCMC : public MetropolisHastings
{
	private:
		std::vector<PsiPrior*> proposaldistributions;
	public:
		DefaultMCMC ( const PsiPsychometric * Model,                                      ///< psychometric function model to sample from
				const PsiData * Data,                                                     ///< data to base inference on
				PsiRandom* proposal                                                       ///< IGNORED
				);
		~DefaultMCMC ( void );
		double acceptance_probability (
				const std::vector<double> &current_theta,
				const std::vector<double> &new_theta );
		void proposePoint ( std::vector<double> &current_theta,
				std::vector<double> &stepwidths,
				PsiRandom * proposal,
				std::vector<double> &new_theta);                                          ///< propose a new sample and save it in new_theta
        void set_proposal(unsigned int i, PsiPrior* proposal){
            delete proposaldistributions.at(i);
            proposaldistributions.at(i) = proposal->clone();
        }
};

class HybridMCMC : public PsiSampler
{
	private:
		PsiRandom* proposal;
		std::vector<double> currenttheta;
		std::vector<double> newtheta;
		std::vector<double> momentum;
		double currentH;
		double newH;
		double energy;
		double newenergy;
		std::vector<double> gradient;
		std::vector<double> currentgradient;
		std::vector<double> stepsizes;
		int Nleapfrog;
		int Naccepted;
		void leapfrog ( void );
	public:
		HybridMCMC (
			const PsiPsychometric * Model,                                                  ///< psychometric function model to sample from
			const PsiData * Data,                                                          ///< data to base inference on
			int Nleap                                                                     ///< number of leapfrog steps to be performed for each sample
			);                                                             ///< initialize the sampler
		~HybridMCMC ( void ) { delete proposal; }
		std::vector<double> draw ( void );                                                ///< draw a sample from the posterior
		void setTheta ( const std::vector<double>& prm );                                 ///< set the current state of the sampler
		std::vector<double> getTheta ( void ) { return currenttheta; }                    ///< get the current state of the sampler
		void setStepSize ( double size, unsigned int param );                             ///< set stepsize of the leapfrog integration for parameter param
		void setStepSize ( const std::vector<double>& sizes );                            ///< set all stepsizes of leapfrog integration for all parameters at once
		double getDeviance ( void );                                                      ///< get the current deviance
		MCMCList sample ( unsigned int N );                                              ///< draw N samples from the posterior
};

/**
 * Model evidence (or marginal likelihood) is given by the following integral
 *
 * \int P(D|theta) * P(theta) d theta.
 *
 * Model evidence can be used to compare models. If E_1 is the evidence for model 1
 * and E_2 is the evidence for model 2. The Bayes Factor for comparison of these two
 * models is
 *
 * BF = E_1/E_2.
 *
 * This can be interpreted as "How much more likeli is model 1 than model 2 given the
 * data.
 *
 * To estimate model evidence, we note that model evidence can be interpreted as the
 * expectation of likelihood under the prior distribution. Thus, this function draws
 * samples from the prior distribution and approximates the above integral as follows:
 *
 * \int P(D|theta) * P(theta) d theta ~~ 1/n \sum_{i=1}^n P(D|theta_i).
 *
 */
double ModelEvidence ( const PsiPsychometric* pmf, const PsiData* data );

/**
 * Bayesian Outlier detection
 *
 * Similar to the procedure proposed by Wichmann & Hill (2001a), outliers are detected
 * by performing a series of model comparisons. The fitted psychometric function model
 * M is compared to a series of models M0,M1,...,Mi, that are modified by incorporating
 * an additional parameter that perfectly fits the i-th block. That is the models Mi
 * have on parameter more than model M (which makes them worse) but they also explain
 * the i-th block perfectly (which makes them better). If the i-th block is not from the
 * same psychometric function as the other blocks, the improvement of adding an additional
 * parameter for the i-th block should be large. That means, model Mi should be preferred
 * to model M. If however, it is plausible that the i-th block comes from the fitted
 * psychometric (defined by all remaining points together), the gain in descriptive power
 * from the additional parameter is weak and model M should be preferred. Such model
 * comparisons can easily be performed using the so called "Bayes Factor", that is
 * the ratio of the evidences of the two models.
 * This function determines a Bayes Factor for each block for the comparison between
 * the model that was based on all blocks and the model that uses a separate parameter
 * for the block of interest. Bayes Factors > 1 favor the global model, Bayes Factors < 1
 * favor the model with a special parameter for block i, indicating that block i is
 * an outlier.
 */
std::vector<double> OutlierDetection ( const PsiPsychometric* pmf, OutlierModel* outl, const PsiData* data );

#endif
