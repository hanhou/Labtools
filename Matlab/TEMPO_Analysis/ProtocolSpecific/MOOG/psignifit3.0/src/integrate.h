#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "psychometric.h"
#include "data.h"
#include "getstart.h"
#include "bootstrap.h"
#include <vector>
#include <algorithm>

/** Produce a linearly spaced grid from xmin to xmax */
std::vector<double> lingrid ( double xmin, double xmax, unsigned int gridsize );

// Numerical integration
void normalize_probability (
		const std::vector<double>& x, ///< nodes at which the density has been evaluated
		std::vector<double>& fx       ///< (unnormalized) density at the respective nodes.
		);  ///< normalize a probability density (in place!)

double numerical_mean (
		const std::vector<double>& x,  ///< nodes at which the density has been evaluated
		const std::vector<double>& fx  ///< normalized density at the respective nodes
		);  ///< determine the 1d mean by evaluating the respective integral numerically

double numerical_variance (
		const std::vector<double>& x,  ///< nodes at which the density has been evaluated
		const std::vector<double>& fx, ///< normalized density at the respective nodes
		double m                       ///< numerically evaluated mean
		);  ///< determine the 1d variance by evaluating the respective integral numerically

// Moment matching
std::vector<double> match_gauss (
		const std::vector<double>& x,  ///< nodes at which the density has been evaluated
		const std::vector<double>& fx  ///< normalized density at the respective nodes
		);  ///< determine parameters of a gaussian that matches with the first and second moments of the sampled density
std::vector<double> match_gamma (
		const std::vector<double>& x,  ///< nodes at which the density has been evaluated
		const std::vector<double>& fx  ///< normalized density at the respective nodes
		);  ///< determine parameters of a gamma distribution that matches with the first and second moments of the sampled density
std::vector<double> match_beta (
		const std::vector<double>& x,  ///< nodes at which the density has been evaluated
		const std::vector<double>& fx  ///< normalized density at the respective nodes
		);  ///< determine parameters of a beta distribution that matches the first and second moments of the sampled density


/** Independent approximation to the posterior distribution */
class PsiIndependentPosterior {
	private:
		unsigned int nparams;
		std::vector<PsiPrior*> fitted_posteriors;
		std::vector< std::vector<double> > grids;
		std::vector< std::vector<double> > margins;
	public:
		PsiIndependentPosterior (
				unsigned int nprm,                    ///< number of parameters in the model
				std::vector<PsiPrior*> posteriors,    ///< determined best fitting distributions
				std::vector< std::vector<double> > x, ///< grids for all parameters
				std::vector< std::vector<double> > fx ///< densities for all parameters
				);
		PsiIndependentPosterior ( const PsiIndependentPosterior& post ) :
			nparams ( post.nparams ), fitted_posteriors ( post.nparams ), grids ( post.grids ), margins ( post.margins )
			{ unsigned int i; for ( i=0; i<nparams; i++ ) fitted_posteriors[i] = post.fitted_posteriors[i]->clone(); }
		~PsiIndependentPosterior ( void ) { unsigned int i; for ( i=0; i<nparams; i++ ) delete fitted_posteriors[i]; }
		PsiPrior *get_posterior ( unsigned int parameter ) { return fitted_posteriors[parameter]->clone(); } ///< get the fitted posterior distribution for a single parameter
		std::vector<double> get_grid ( unsigned int parameter ) { return grids[parameter]; } ///< get the grid points of a single parameter
		std::vector<double> get_margin ( unsigned int parameter ) { return margins[parameter]; } ///< get the (marginal) density of a single parameter
};

PsiIndependentPosterior independent_marginals (
		const PsiPsychometric *pmf,    ///< psychometric function model
		const PsiData *data            ///< dataset
		);  ///< determine an approximation to the posterior distribution that approximates the posterior as a product of independent distributions for all parameters

MCMCList sample_posterior (
		const PsiPsychometric *pmf,    ///< psychometric function model
		const PsiData *data,           ///< dataset
		PsiIndependentPosterior& post, ///< posterior approximation
		unsigned int nsamples=600,     ///< number of samples to be drawn
		unsigned int propose=25        ///< oversampling factor for the proposals
		);   ///< sample from the posterior using sampling importance resampling

void sample_diagnostics (
		const PsiPsychometric *pmf,   ///< psychometric function model
		const PsiData *data,          ///< dataset
		MCMCList *samples             ///< parameter samples
		);  ///< calculate sample diagnostics

#endif
