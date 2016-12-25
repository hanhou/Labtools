/**
 * psignifit-mapestimate:  mapestimation of parameters for psychometric functions
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

int main ( int argc, char ** argv ) {
	// Analyze the command line
	cli_parser parser ( "psignifit-mapestimate [options] <file> [ <file> ... ]" );
	parser.add_option ( "-c", "psignifit core object to be used", "mw0.1" );
	parser.add_option ( "-s", "psignifit sigmoid object to be used", "logistic" );
	parser.add_option ( "-prior1", "prior for the first parameter (alpha,a,m,...)", "None" );
	parser.add_option ( "-prior2", "prior for the second parameter (beta,b,w,...)", "None" );
	parser.add_option ( "-prior3", "prior for the third parameter (lambda)", "Uniform(0,.1)" );
	parser.add_option ( "-prior4", "prior for the fourth parameter (gamma)", "Uniform(0,.1)" );
	parser.add_option ( "-nafc",   "number of response alternatives in forced choice designs (set this to 1 for yes-no tasks)", "2" );
	parser.add_option ( "-o",      "write output to this file", "stdout" );
	parser.add_option ( "-cuts",   "cuts to be determined", "0.25,0.50,0.75" );
	parser.add_switch ( "-v", "display status messages", false );
	parser.add_switch ( "-e", "In yes-no tasks: set gamma==lambda", false );
	parser.add_switch ( "--matlab", "format output to be parsable by matlab", false );

	parser.parse_args ( argc, argv );

	// Set up the environment
	bool verbose ( parser.getOptSet ( "-v" ) ), pmfshown ( false );
	PsiData *data;
	PsiPsychometric *pmf;
	PsiOptimizer * opt;
	std::vector<double> theta;
	std::vector<double> cuts (getCuts ( parser.getOptArg("-cuts") ) );
	unsigned int i;

	// Where do we want output to go
	FILE * ofile;
	if ( !(parser.getOptArg ( "-o" ).compare( "stdout" )) ) {
		if ( verbose ) std::cerr << "Writing results to stdout\n";
		ofile = stdout;
	} else 
		ofile = fopen ( parser.getOptArg ( "-o" ).c_str(), "w" );

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
	}

	std::string fname;
	fname = parser.popArg ();
	if ( fname == "" ) {
		std::cerr << "No input file given --- aborting!\n";
		exit ( -1 );
	}

	while ( fname != "" ) {
		if ( verbose ) std::cerr << "Analyzing input file '" << fname << "'\n   ";

		data = allocateDataFromFile ( fname, atoi ( parser.getOptArg ( "-nafc" ).c_str() ) );

		if ( verbose ) std::cerr << "Read " << data->getNblocks() << " blocks ";

		// Here, we set up the psychometric function model
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

		// Perform the optimization
		opt = new PsiOptimizer ( pmf, data );
		theta = opt->optimize ( pmf, data );
		// Replace cuts with the derived thresholds at the cuts
		for ( i=0; i<cuts.size(); i++ ) {
			cuts[i] = pmf->getThres ( theta, cuts[i] );
		}

		if ( verbose ) std::cerr << "\n";

		// Print output
		print ( theta,                          parser.getOptSet ( "--matlab" ), "params_estimate", ofile );
		print_fisher ( pmf, theta, data, ofile, parser.getOptSet ( "--matlab" ) );
		print ( cuts,                           parser.getOptSet ( "--matlab" ), "thres", ofile );
		print ( pmf->deviance ( theta, data ),  parser.getOptSet ( "--matlab" ), "deviance", ofile );

		fname = parser.popArg();

		// Clean up
		delete data;
		delete pmf;
		delete opt;
	}

	return 0;
}
