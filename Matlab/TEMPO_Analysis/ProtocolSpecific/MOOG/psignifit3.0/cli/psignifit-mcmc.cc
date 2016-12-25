/**
 * psignifit-mcmc: bayesian inference for psychometric functions
 *
 * This file is part of the command line interface to psignifit.
 *
 * See COPYING file distributed along with the psignifit package for the copyright and
 * license terms
 */
#include "../src/psipp.h"
#include "cli.h"

#include "cli_utilities.h"

#include <cstdio>
#include <algorithm>

double prctile ( std::vector<double> x, double p ) {
	sort ( x.begin(), x.end() );
	unsigned int pos ( floor(p*x.size()) );
	return x[pos];
}

int main ( int argc, char ** argv ) {
	// Parse command line
	cli_parser parser ( "psignifit-mcmc [options] <file> [ <file> ... ]" );
	parser.add_option ( "-c",           "psignifit core object to be used", "mw0.1" );
	parser.add_option ( "-s",           "psignifit sigmoid object to be used", "logistic" );
	parser.add_option ( "-prior1",      "prior for the first parameter (alpha,a,m,...)", "None" );
	parser.add_option ( "-prior2",      "prior for the second parameter (beta,b,w,...)", "None" );
	parser.add_option ( "-prior3",      "prior for the third parameter (lambda)", "Uniform(0,.1)" );
	parser.add_option ( "-prior4",      "prior for the fourth parameter (gamma)", "Uniform(0,.1)" );
	parser.add_option ( "-nafc",        "number of response alternatives in forced choice designs (set this to 1 for yes-no tasks)", "2" );
	parser.add_option ( "-nsamples",    "number of markov chain monte carlo samples to be generated","2000" );
	parser.add_option ( "-o",           "write output to this file", "stdout" );
	parser.add_option ( "-cuts",        "cuts to be determined", "0.25,0.50,0.75" );
	parser.add_option ( "-proposal",    "standard deviations of the proposal distribution (or name of file with pilot samples)", "0.1,0.1,0.01" );
	parser.add_option ( "-start",       "starting values for the sampling process", "mapestimate" );
	parser.add_switch ( "-v",           "display status messages", false );
	parser.add_switch ( "--summary",    "write a short summary to stdout" );
	parser.add_switch ( "-e",           "In yes-no tasks: set gamma==lambda", false );
	parser.add_switch ( "-generic",     "Use generic metropolis instead of the default standard metropolis hastings", false );
	parser.add_switch ( "--matlab",     "format output to be parsable by matlab", false );

	parser.parse_args ( argc, argv );

	// Set up the most important data
	bool verbose ( parser.getOptSet ( "-v" ) ), pmfshown ( false ), summary ( parser.getOptSet( "--summary" ) ), generic ( parser.getOptSet ( "-generic" ) );
	unsigned int i,j, ncuts, nparams, nblocks;
	PsiData                    *data;
	PsiPsychometric            *pmf;
	PsiOptimizer               *opt;
	PsiSampler                 *sampler;
	std::vector<double>         theta;
	std::vector<double>         cuts (getCuts ( parser.getOptArg("-cuts") ) );
	                            ncuts = cuts.size ();
	PsiMClist                  *pilotsample = NULL;
	char                        sline[80];
	MCMCList                   *mcmc_list;
	unsigned int                nsamples ( atoi ( parser.getOptArg("-nsamples").c_str() ) );
	double                      th;
	double 						sl;
	double 						th_m;
	double 						sl_m;	
	double                      meanestimate;
	double                      bayesian_p;
	GaussRandom                 proposal;
	setSeed (0);

	// We might either want to read the stepwidths or a pilot sample
	std::vector<double>         stepwidths;
	std::fstream                pilotfile ( parser.getOptArg ( "-proposal" ).c_str() );
	if ( pilotfile.fail() )
		stepwidths = getCuts ( parser.getOptArg("-proposal") );
	else {
		// read data from file
		if ( verbose ) std::cerr << "Reading pilot data...";
		while ( strcmp(sline,"# mcestimates") ) {
			pilotfile.getline ( sline, 80 );
		}
		j = pilotfile.tellg();
		nblocks = 1;
		pilotfile.getline ( sline, 80 );
		nparams = 0;
		i = std::string ( sline ).find('.');
		while ( i!=std::string::npos ) {
			nparams ++;
			i = std::string ( sline ).find ('.',i+1);
		}
		if ( verbose ) std::cerr << "nparams: " << nparams;
		while ( strcmp(sline,"") ) {
			nblocks++;
			pilotfile.getline ( sline, 80 );
		}
		pilotfile.clear(); // If we passed the end of the file
		pilotfile.seekg ( j, std::ios::beg );
		theta = std::vector<double> ( nparams );
		pilotsample = new PsiMClist ( nblocks, nparams );
		for ( i=0; i<nblocks; i++ ) {
			for ( j=0; j<nparams; j++ ) {
				pilotfile >> theta[j];
			}
			pilotsample->setEst ( i, theta, 0 );
		}
		// In case we use MH-Sampling: determine stepwidths as averages
		if ( verbose ) std::cerr << " sd: ";
		stepwidths = std::vector<double> ( nparams );
		for ( j=0 ; j<nparams; j++ ) {
			stepwidths[j] = 0;
			theta[j] = pilotsample->getMean ( j );
			for ( i=0; i<nblocks; i++ ) {
				stepwidths[j] += (pilotsample->getEst ( i, j )-theta[j])*(pilotsample->getEst ( i, j )-theta[j]);
			}
			stepwidths[j] /= pilotsample->getNsamples()-1;
			stepwidths[j] = sqrt(stepwidths[j]);
			if ( verbose ) std::cerr << stepwidths[j] << " ";
		}
		if ( verbose ) std::cerr << "\n";
	}

	// Contents of the mcmc lists
	std::vector< std::vector<double> >  mcthres ( nsamples, cuts );
	std::vector< std::vector<double> >  mcslopes ( nsamples, cuts );
	std::vector< std::vector<double> > *mcestimates;
	std::vector< std::vector<int> >    *mcdata;
	std::vector<double>                 mcdeviance ( nsamples );
	std::vector<double>                 mcRpd      ( nsamples );
	std::vector<double>                 mcRkd      ( nsamples );
	std::vector<double>                 ppdeviance ( nsamples );
	std::vector<double>                 ppRpd      ( nsamples );
	std::vector<double>                 ppRkd      ( nsamples );
	std::vector<double>                 dummydata  ( nsamples/2 );
	std::vector<double>   thresholds ( ncuts );
	std::vector<double>   slopes     ( ncuts );
	std::vector<double>                *influential;
	std::vector<int>                   *outliers;
	std::vector<double>                *ci_lower;
	std::vector<double>                *ci_upper;
	std::vector<double>                *devianceresiduals;

	// Use matlabformat?
	bool matlabformat ( parser.getOptSet ( "--matlab" ) );

	// Get the output file
	FILE * ofile;
	if ( !(parser.getOptArg ( "-o" ).compare( "stdout" )) ) {
		if ( verbose ) std::cerr << "Writing results to stdout\n";
		ofile = stdout;
	} else 
		ofile = fopen ( parser.getOptArg ( "-o" ).c_str(), "w" );

	// Write some status messages
	if (verbose) {
		std::cerr << "core:    " << parser.getOptArg ( "-c" ) << "\n";
		std::cerr << "sigmoid: " << parser.getOptArg ( "-s" ) << "\n";
		std::cerr << "cuts:    ";
		for (i=0; i<cuts.size(); i++) std::cerr << cuts[i] << " ";
		std::cerr << "\n";
		std::cerr << "priors:\n";
		std::cerr << "   prm1: " << parser.getOptArg ( "-prior1" ) << "\n";
		std::cerr << "   prm2: " << parser.getOptArg ( "-prior2" ) << "\n";
		std::cerr << "   prm3: " << parser.getOptArg ( "-prior3" ) << "\n";
		if ( atoi (parser.getOptArg("-nafc").c_str()) < 2 ) std::cerr << "   prm4: " << parser.getOptArg ( "-prior4" ) << "\n";
		std::cerr << "generic mcmc: " << (parser.getOptSet("-generic")?"yes":"no") << "\n";
		std::cerr << "number of mcmc samples: " << nsamples << "\n";
		if ( parser.getOptSet ( "-e" ) )
			std::cerr << "gamma==lambda\n";
		std::cerr << "stepwidths:\n";
		std::cerr << "   s1: " << stepwidths[0] << "\n";
		std::cerr << "   s2: " << stepwidths[1] << "\n";
		std::cerr << "   s3: " << stepwidths[2] << "\n";
		if ( stepwidths.size()==4 ) std::cerr << "   s4: " << stepwidths[3] << "\n";
	}

	std::string fname;
	fname = parser.popArg ();
	if ( fname == "" ) {
		std::cerr << "No input file given --- aborting!\n";
		exit ( -1 );
	}

	while ( fname != "" ) {
		if ( verbose ) std::cerr << "Analyzing input file '" << fname << "'\n   ";

		// Get the data
		data = allocateDataFromFile ( fname, atoi ( parser.getOptArg ( "-nafc" ).c_str() ) );
		nblocks = data->getNblocks();

		if ( verbose ) std::cerr << "Read " << nblocks << " blocks ";

		// Get the psychometric function model
		pmf  = allocatePsychometric ( parser.getOptArg ( "-c" ),
			parser.getOptArg ( "-s" ),
			atoi ( parser.getOptArg ( "-nafc" ).c_str() ),
			data,
			verbose && !pmfshown);
		pmfshown = true;
		if ( parser.getOptSet ( "-e" ) ) pmf->setgammatolambda();
		setPriors ( pmf,
				parser.getOptArg ( "-prior1" ),
				parser.getOptArg ( "-prior2" ),
				parser.getOptArg ( "-prior3" ),
				parser.getOptArg ( "-prior4" ) );
		nparams = pmf->getNparams();

		// Determine starting value
		opt = new PsiOptimizer ( pmf, data );
		theta = opt->optimize ( pmf, data );
		
		// Set up the sampler
		if ( generic ) {
			sampler = new GenericMetropolis ( pmf, data, &proposal );
			if ( pilotsample != NULL ) {
				((GenericMetropolis*)sampler)->findOptimalStepwidth ( *pilotsample );
			} else {
				sampler->setStepSize ( stepwidths );
			}
		} else {
			sampler = new MetropolisHastings ( pmf, data, &proposal );
			sampler->setStepSize ( stepwidths );
		}
		if ( parser.getOptArg ( "-start" ) == "mapestimate" ) {
			((MetropolisHastings*)sampler)->setTheta ( theta );
		} else {
			((MetropolisHastings*)sampler)->setTheta ( getCuts ( parser.getOptArg ( "-start" ) ) );
		}

		// Sample
		if ( verbose ) {
			std::cerr << "Starting sampling ...";
			std::cerr.flush();
		}
		mcmc_list = new MCMCList ( sampler->sample ( nsamples ) );

		if ( verbose ) std::cerr << " Done \n";

		// These might change during analysis and have to be allocated for each file
		mcestimates = new std::vector< std::vector<double> > (nsamples);
		mcdata      = new std::vector< std::vector<int> >    (nsamples);
		influential = new std::vector<double>                (nblocks);
		outliers    = new std::vector<int>                   (nblocks);
		ci_lower    = new std::vector<double>                (nparams);
		ci_upper    = new std::vector<double>                (nparams);

		// Now store everything that is related to inteval estimation of parameters
		for ( i=0; i<nsamples; i++ ) {
			(*mcestimates)[i] = mcmc_list->getEst ( i );
			(*mcdata)[i]      = mcmc_list->getppData ( i );
			for ( j=0; j<ncuts; j++ ) {
				mcthres[i][j]  = pmf->getThres ( (*mcestimates)[i], cuts[j] );
				mcslopes[i][j] = pmf->getSlope ( (*mcestimates)[i], cuts[j] );
			}
		}
		for ( i=0; i<nblocks; i++ ) {
			(*influential)[i] = 0;
			for ( j=nsamples/2; j<nsamples; j++ ) {
				(*influential)[i] += mcmc_list->getlogratio (j,i);
			}
		}

		// Write a summary of the parameter estimation if requested.
		if ( summary ) {
			std::cerr << "Parameter estimates:\n";
			std::cerr << "--------------------\n";
			for ( i=0; i<nparams; i++ ) {
				meanestimate = 0;
				for ( j=0; j<nsamples/2; j++ ) {
					dummydata[j] = (*mcestimates)[nsamples/2+j][i];
					meanestimate += dummydata[j];
				}
				meanestimate /= nsamples;
				// Why are the parameter estiamtes so strange?
				std::cerr << "parameter" << i+1 << " = " << meanestimate << "\tCI_95 = (" << prctile(dummydata,.025) << "," << prctile(dummydata,.975) << ")\n";
				theta[i] = meanestimate;
			}
			std::cerr << "\n";
			std::cerr << "Threshold estimates:\n";
			std::cerr << "--------------------\n";
			for ( i=0; i<ncuts; i++ ) {
				th = pmf->getThres ( theta, cuts[i] );

				for ( j=0; j<nsamples/2; j++ ) {
					dummydata[j] = pmf->getThres ( (*mcestimates)[nsamples/2+j], cuts[i] );
				}
				std::cerr << "Threshold(" << cuts[i] << ") = " << th << "\tCI_95 = ("
					<< prctile(dummydata,.025) << ","
					<< prctile(dummydata,.975) << ") ";
				for ( j=0; j<nsamples/2; j++ ) {
					dummydata[j] = pmf->getSlope ( (*mcestimates)[nsamples/2+j], cuts[i] );
				}
				std::cerr << "Slope(" << cuts[i] << ") = " << pmf->getSlope ( theta, th ) << "\tCI_95 = ("
					<< prctile(dummydata,.025) << ","
					<< prctile(dummydata,.975) << ")\n";
			}
		}

		
		for ( i=0; i<ncuts; i++ ) {
		th_m = pmf->getThres ( theta, cuts[i] );
		thresholds[i] = th_m;
		sl_m = pmf->getSlope ( theta, th_m );
		slopes[i] = sl_m;
		}
		
		// Now store everything related to goodness of fit
		for ( i=0; i<nsamples; i++ ) {
			mcdeviance[i] = mcmc_list->getdeviance(i);
			mcRpd[i]      = mcmc_list->getRpd(i);
			mcRkd[i]      = mcmc_list->getRkd(i);
			ppdeviance[i] = mcmc_list->getppDeviance(i);
			ppRpd[i]      = mcmc_list->getppRpd(i);
			ppRkd[i]      = mcmc_list->getppRkd(i);
		}

		// Write a summary of the goodness of fit statistics if requested
		if ( summary ) {
			devianceresiduals = new std::vector<double> ( pmf->getDevianceResiduals ( theta, data ) );
			std::cerr << "\n";
			std::cerr << "Goodness of fit statistics:\n";
			std::cerr << "---------------------------\n";
			bayesian_p = 0; for ( i=nsamples/2; i<nsamples; i++ ) bayesian_p += ppdeviance[i] > mcdeviance[i]; bayesian_p /= nsamples/2;
			std::cerr << "Deviance: " << pmf->deviance ( theta, data ) <<             "\tbayesian_p: " << bayesian_p << "\n";
			bayesian_p = 0; for ( i=nsamples/2; i<nsamples; i++ ) bayesian_p += ppRpd[i] > mcRpd[i]; bayesian_p /= nsamples/2;
			std::cerr << "Rpd:      " << pmf->getRpd (
					*devianceresiduals, theta, data ) << "\tbayesian_p: " << bayesian_p << ")\n";
			bayesian_p = 0; for ( i=nsamples/2; i<nsamples; i++ ) bayesian_p += ppRkd[i] > mcRkd[i]; bayesian_p /= nsamples/2;
			std::cerr << "Rkd:      " << pmf->getRkd (
					*devianceresiduals, data ) <<        "\tbayesian_p: " << bayesian_p << ")\n";
			delete devianceresiduals;
		}

		// Now store everything in the output file.
		print ( *mcdata,      matlabformat, "mcdata",      ofile );
		print ( *mcestimates, matlabformat, "mcestimates", ofile );
		print ( mcdeviance,   matlabformat, "mcdeviance",  ofile );
		print ( ppdeviance,   matlabformat, "ppdeviance",  ofile );
		print ( mcthres,      matlabformat, "mcthres",     ofile );
		print ( mcslopes,     matlabformat, "mcslopes",    ofile );
		print ( mcRpd,        matlabformat, "mcRpd",       ofile );
		print ( ppRpd,        matlabformat, "ppRpd",       ofile );
		print ( mcRkd,        matlabformat, "mcRkd",       ofile );
		print ( ppRkd,        matlabformat, "ppRkd",       ofile );
		print ( *influential, matlabformat, "influential", ofile );
		print ( thresholds,   matlabformat, "thresholds",  ofile );
		print ( slopes,       matlabformat, "slopes",      ofile );

		// Get the next input file (if there is one)
		fname = parser.popArg();

		// Clean up
		delete mcestimates;
		delete mcdata;
		delete influential;
		delete outliers;
		delete data;
		delete pmf;
		delete opt;
	}

	if (pilotsample!=NULL) delete pilotsample;

	return 0;
}

