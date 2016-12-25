#include "psipp.h"
#include <R.h>
#include <Rmath.h>
#include <vector>
#include <cstdio>

#include <iostream>

PsiSigmoid * determine_sigmoid ( char *sigmoid ) {
	if ( !strcmp(sigmoid,"logistic") ) {
		return new PsiLogistic ();
	} else if ( !strcmp(sigmoid,"exp") ) {
		return new PsiExponential ();
	} else if ( !strcmp(sigmoid,"gauss") ) {
		return new PsiGauss ();
	} else if ( !strcmp(sigmoid,"gumbel_l") || !strcmp(sigmoid,"gumbel") || !strcmp(sigmoid,"lgumbel") ) {
		return new PsiGumbelL ();
	} else if ( !strcmp(sigmoid,"gumbel_r") || !strcmp(sigmoid,"rgumbel") ) {
		return new PsiGumbelL();
	} else if ( !strcmp(sigmoid,"cauchy") ) {
		return new PsiCauchy();
	} else {
		Rprintf ( "WARNING: no valid sigmoid!" );
		return NULL;
	}
}

PsiCore * determine_core ( char *core, const PsiSigmoid* Sigmoid, const PsiData* data ) {
	double dummy;
	if ( !strcmp(core,"ab") ) {
		return  new abCore();
	} else if ( !strncmp(core,"mw",2) ) {
		sscanf(core, "mw%lf", &dummy);
		return  new mwCore(Sigmoid->getcode(), dummy);
	} else if ( !strcmp(core,"linear") ) {
		return  new linearCore();
	} else if ( !strcmp(core,"log") ) {
		return  new logCore(data);
	} else if ( !strcmp(core,"poly") ) {
		return  new polyCore(data);
	} else if ( !strcmp(core,"weibull") ) {
		return  new weibullCore(data);
	} else {
		Rprintf ( "WARNING: no valid core!" );
		return NULL;
	}
}

PsiPrior * determine_prior ( char *prior ) {
	double prm1,prm2;
	if ( !strncmp(prior,"Beta",4) ) {
		sscanf(prior,"Beta(%lf,%lf)",&prm1,&prm2);
		return new BetaPrior(prm1,prm2);
	} else if ( !strncmp(prior,"Gamma",5) ) {
		sscanf(prior,"Gamma(%lf,%lf)",&prm1,&prm2);
		return new GammaPrior(prm1,prm2);
	} else if ( !strncmp(prior,"Gauss",5) ) {
		sscanf(prior,"Gauss(%lf,%lf)",&prm1,&prm2);
		return new GaussPrior(prm1,prm2);
	} else if ( !strncmp(prior,"Uniform",7) ) {
		sscanf(prior,"Uniform(%lf,%lf)",&prm1,&prm2);
		return new UniformPrior(prm1,prm2);
	} else {
		return new PsiPrior();
	}
}

PsiData * determine_data ( double *x, int *k, int *n, int *K, int *nafc ) {
	std::vector<double> stimulus_intensities ( *K );
	std::vector<int>    number_of_trials     ( *K );
	std::vector<int>    number_of_correct    ( *K );
	int i;

	for ( i=0; i<*K; i++ ) {
		stimulus_intensities[i] = x[i];
		number_of_trials[i]     = n[i];
		number_of_correct[i]    = k[i];
	}

	return new PsiData ( stimulus_intensities, number_of_trials, number_of_correct, *nafc );
}

void get_fitting_setup ( double *x, int *k, int *n, int *K,
		char **sigmoid, char **core, int *nafc, int *nparams, char **priors,
		PsiData** dataout, PsiPsychometric** pmfout )
{
	int i;
	*dataout = determine_data ( x, k, n, K, nafc );

	PsiSigmoid *Sigmoid = determine_sigmoid ( *sigmoid );
	if (Sigmoid==NULL) { delete *dataout; throw -1; }

	PsiCore *Core = determine_core ( *core, Sigmoid, *dataout );
	if (Core==NULL) { delete *dataout; delete Sigmoid; throw -1; }

	*pmfout = new PsiPsychometric ( *nafc, Core, Sigmoid );
	if (*nparams != (*pmfout)->getNparams() ) {
		Rprintf ( "WARNING: output vector length does not match number of parameters!" );
		delete dataout;
		delete Sigmoid;
		delete Core;
		delete *pmfout;
		throw -1;
	}

	for ( i=0; i<*nparams; i++ ) {
		(*pmfout)->setPrior ( i, determine_prior ( priors[i] ) );
	}
}

extern "C" {
////////////////////////////////////
// Functions go here
////////////////////////////////////

void mapestimate (
		double *x,        // stimulus intensities
		int *k,           // response counts of correct- (nAFC) or Yes-responses (Yes/No)
		int *n,           // numbers of trials per block
		int *K,           // number of blocks
		char **sigmoid,   // the sigmoid to be used
		char **core,      // core description
		int *nafc,        // number of alternatives in the task (a value < 2 indicates Yes/No)
		char **priors,    // priors
		int *nparams,     // number of parameters
		double *estimate, // output array for the estimated values
		double *deviance, // output: deviance
		double *Rpd,      // output: correlation between model prediction and deviance residuals
		double *Rkd       // output: correlation between block index and deviance residuals
		) {
	int i;
	PsiData * data;
	PsiPsychometric *pmf;

	try {
		get_fitting_setup ( x, k, n, K, sigmoid, core, nafc, nparams, priors, &data, &pmf );
	} catch (int) {
		return;
	}

	std::vector<double> startest = pmf->getStart ( data );
	PsiOptimizer *opt = new PsiOptimizer ( pmf, data );
	std::vector<double> est = opt->optimize ( pmf, data );
	delete opt;


	for ( i=0; i<*nparams; i++ )
		estimate[i] = est[i];

	std::vector<double> devianceresiduals ( pmf->getDevianceResiduals ( est, data ) );
	*deviance = pmf->deviance(est,data);
	*Rpd      = pmf->getRpd ( devianceresiduals, est, data );
	*Rkd      = pmf->getRkd ( devianceresiduals, data );

	delete data;
	delete pmf;

	return;
}

void performbootstrap (
		double *x,          // stimulus intensities
		int *k,             // response counts of correct- (nAFC) or Yes-responses (Yes/No)
		int *n,             // numbers of trials per block
		int *K,             // number of blocks
		char **sigmoid,     // the sigmoid to be used
		char **core,        // core description
		int *nafc,          // number of alternatives in the task (a value < 2 indicates Yes/No)
		char **priors,      // priors
		double *generating, // generating probabilites (if generating[0]==-999, non parametric bootstrap is performed)
		int *nparams,       // number of parameters
		int *nsamples,      // number of bootstrap samples to be drawn
		double *cuts,       // cuts at which the thresholds should be determined
		int *ncuts,         // number of cuts
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		int *bdata,         // output: bootstrap samples (needs to be 'reshaped')
		double *bestimates, // output: bootstrap estimates (needs to be reshaped)
		double *bdeviances, // output: deviances of the bootstrap samples
		double *bRpd,       // output: correlations between model prediction and deviance residuals for all bootstrap samples
		double *bRkd,       // output: correlations between block index and deviance residuals for all bootstrap samples
		double *bthres,     // output: thresholds of the bootstrap samples
		double *influential,// output: influence of the different blocks
		double *acc,        // output: acceleration constants
		double *bias,       // output: bias correction constants
		double *thresholdci // output: 95% confidence intervals for thresholds
		) {
	int i,j;
	PsiData * data;
	PsiPsychometric * pmf;
	std::vector<double> *start;
	bool parametric;

	try {
		get_fitting_setup ( x, k, n, K, sigmoid, core, nafc, nparams, priors, &data, &pmf );
	} catch (int) {
		return;
	}

	std::vector<double> Cuts ( *ncuts );
	for ( i=0; i<*ncuts; i++ ) Cuts[i] = cuts[i];

	if (generating[0]==-999) {
		start = NULL;
		parametric = false;
	} else {
		start = new std::vector<double> (*nparams);
		for ( i=0; i<*nparams; i++ ) {
			(*start)[i] = generating[i];
		}
		parametric = true;
	}

	BootstrapList bslist ( bootstrap (*nsamples,data,pmf,Cuts,start,true, parametric) );
	JackKnifeList jack   ( jackknifedata ( data, pmf ) );

	for ( i=0; i<*nsamples; i++ ) {
		for ( j=0; j<*K; j++ ) {
			bdata[*K*i + j] = bslist.getData(i)[j];
		}
		for ( j=0; j<*nparams; j++ ) {
			bestimates[*nparams*i + j] = bslist.getEst(i, j);
		}
		bdeviances[i] = bslist.getdeviance(i);
		bRpd[i]       = bslist.getRpd(i);
		bRkd[i]       = bslist.getRkd(i);
		for ( j=0; j<*ncuts; j++ ) {
			bthres[*ncuts*i + j] = bslist.getThres_byPos (i,j);
		}
	}
	for ( j=0; j<*ncuts; j++ ) {
		acc[j] = bslist.getAcc(j);
		bias[j] = bslist.getBias(j);
		thresholdci[3*j]   = bslist.getThres(0.025, j);
		thresholdci[3*j+1] = bslist.getThres(0.5, j);
		thresholdci[3*j+2] = bslist.getThres(0.975, j);
	}

	std::vector<double> ci_lower ( *nparams );
	std::vector<double> ci_upper ( *nparams );
	for ( j=0; j<*nparams; j++ ) {
		ci_lower[j] = bslist.getPercentile(0.025,j);
		ci_upper[j] = bslist.getPercentile(0.975,j);
	}
	for ( j=0; j<*K; j++ ) {
		influential[j] = jack.influential ( j, ci_lower, ci_upper );
	}

	delete data;
	delete pmf;
	delete start;

	return;
}

void performmcmc (
		double *x,          // stimulus intensities
		int *k,             // response counts of correct- (nAFC) or Yes-responses (Yes/No)
		int *n,             // numbers of trials per block
		int *K,             // number of blocks
		char **sigmoid,     // the sigmoid to be used
		char **core,        // core description
		int *nafc,          // number of alternatives in the task (a value < 2 indicates Yes/No)
		char **priors,      // priors
		double *proposal,   // standard deviances of the proposal distributions
		double *start,      // starting value of the sampler
		int *nparams,       // number of parameters
		int *nsamples,      // number of bootstrap samples to be drawn
		double *cuts,       // cuts at which the thresholds should be determined
		int *ncuts,         // number of cuts
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		double *mcmcestimates, // output: bootstrap estimates (needs to be reshaped)
		double *mcmcdeviances, // output: deviances of the bootstrap samples
		double *mcmcRpd,       // output: correlations between model prediction and deviance residuals for mcmc sammples
		double *mcmcRkd,       // output: correlations between block index and deviance residuals for all mcmc samples
		int *ppdata,           // output: bootstrap samples (needs to be 'reshaped')
		double *ppRpd,         // output: correlations between model prediction and deviance residuals for all bootstrap samples
		double *ppRkd,         // output: correlations between block index and deviance residuals for all bootstrap samples
		double *ppdeviances,   // output: acceleration constants
		double *logpratio      // output: logposterior ratio (for determination of influential observations)
		) {
	int i,j;
	PsiData * data;
	PsiPsychometric * pmf;

	try {
		get_fitting_setup ( x, k, n, K, sigmoid, core, nafc, nparams, priors, &data, &pmf );
	} catch (int) {
		return;
	}

	MetropolisHastings S ( pmf, data, new GaussRandom () );
	for ( i=0; i<*nparams; i++ ) S.setstepsize ( proposal[i], i );

	MCMCList mcmclist ( S.sample( *nsamples ) );

	for ( i=0; i<*nsamples; i++ ) {
		for ( j=0; j<*nparams; j++ ) {
			mcmcestimates[*nparams*i + j] = mcmclist.getEst (i, j);
		}
		mcmcdeviances[i] = mcmclist.getdeviance ( i );
		mcmcRpd[i] = mcmclist.getRpd ( i );
		mcmcRkd[i] = mcmclist.getRkd ( i );

		for ( j=0; j<*K; j++ ) {
			ppdata[*K*i + j] = mcmclist.getppData ( i, j );
		}
		ppRpd[i] = mcmclist.getppRpd ( i );
		ppRkd[i] = mcmclist.getppRkd ( i );
		ppdeviances[i] = mcmclist.getppDeviance ( i );

		for ( j=0; j<*K; j++ ) {
			logpratio[*K*i + j] = mcmclist.getlogratio ( i, j );
		}
	}
}

void getdiagnostics (
		double *x,        // stimulus intensities
		int *k,           // response counts of correct- (nAFC) or Yes-responses (Yes/No)
		int *n,           // numbers of trials per block
		int *K,           // number of blocks
		char **sigmoid,   // the sigmoid to be used
		char **core,      // core description
		int *nafc,        // number of alternatives in the task (a value < 2 indicates Yes/No)
		char **priors,    // priors
		int *nparams,     // number of parameters
		double *cuts,     // cuts at which the thresholds should be determined
		int *ncuts,       // number of cuts
		double *estimate, // array of the estimated values
		double *deviance, // output: deviance
		double *Rpd,      // output: correlation between model prediction and deviance residuals
		double *Rkd,      // output: correlation between block index and deviance residuals
		double *thres,    // output: thresholds at the respective cuts
		double *devianceresiduals // output: deviance residuals
		) {
	int i;
	PsiData * data;
	PsiPsychometric * pmf;
	std::vector<double> theta ( *nparams );
	for ( i=0; i<*nparams; i++ ) theta[i] = estimate[i];

	try {
		get_fitting_setup ( x, k, n, K, sigmoid, core, nafc, nparams, priors, &data, &pmf );
	} catch (int) {
		return;
	}

	std::vector<double> dr ( pmf->getDevianceResiduals ( theta, data ) );

	*deviance = pmf->deviance ( theta, data );
	*Rpd      = pmf->getRpd   ( dr, theta, data );
	*Rkd      = pmf->getRkd   ( dr, data );

	for ( i=0; i<*K; i++ ) {
		devianceresiduals[i] = dr[i];
	}

	for ( i=0; i<*ncuts; i++ ) {
		thres[i] = pmf->getThres ( theta, cuts[i] );
	}

	delete data;
	delete pmf;

	return;
}

void pmfevaluate (
		double *x,        // x values at which the psychometric function should be evaluated
		int *lenx,        // length of x
		double *params,   // parameter values
		int *nparams,  // number of parameters
		char **sigmoid,   // the sigmoid to be used
		char **core,      // core description
		int *nafc,        // number of alternatives in the task (a value < 2 indicates Yes/No)
		double *Fx        // output: output of the psychometric function
		) {
	int i;
	std::vector<double> xx ( 2 );
	std::vector<int> kk ( 2 );
	std::vector<int> nn ( 2 );
	PsiData * dummydata = new PsiData ( xx, kk, nn, *nafc );
	PsiPsychometric * pmf;
	std::vector<double> theta ( *nparams );
	for ( i=0; i<*nparams; i++ ) theta[i] = params[i];

	PsiSigmoid *Sigmoid = determine_sigmoid ( *sigmoid );
	PsiCore *   Core    = determine_core ( *core, Sigmoid, dummydata );
	delete dummydata;
	pmf = new PsiPsychometric ( *nafc, Core, Sigmoid );

	for ( i=0; i<*lenx; i++  ) {
		Fx[i] = pmf->evaluate ( x[i], theta );
	}

	delete pmf;

	return;
}


// Some more?
}
