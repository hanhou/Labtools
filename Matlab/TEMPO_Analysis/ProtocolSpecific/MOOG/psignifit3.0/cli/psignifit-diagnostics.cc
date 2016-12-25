/**
 * psignifit-diagnostics: diagnostic information for psychometric functions
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
	cli_parser parser ( "psignifit-diagnostics [options] <file> [ <file> ... ]" );
	parser.add_option ( "-c", "psignifit core object to be used", "mw0.1" );
	parser.add_option ( "-s", "psignifit sigmoid object to be used", "logistic" );
	parser.add_option ( "-nafc",   "number of response alternatives in forced choice designs (set this to 1 for yes-no tasks)", "2" );
	parser.add_option ( "-o",      "write output to this file", "stdout" );
	parser.add_option ( "-cuts",   "cuts to be determined", "0.25,0.50,0.75" );
	parser.add_option ( "-params", "parameters to be evaluated", "4.0,2.0,0.02" );
	parser.add_switch ( "-v", "display status messages", false );
	parser.add_switch ( "-e", "In yes-no tasks: set gamma==lambda", false );
	parser.add_switch ( "--matlab", "format output to be parsable by matlab", false );

	parser.parse_args ( argc, argv );

	// Set up the environment
	bool verbose ( parser.getOptSet ( "-v" ) ), pmfshown ( false );
	PsiData *data;
	PsiPsychometric *pmf;
	std::vector<double> theta (getCuts ( parser.getOptArg("-params") ) );
	std::vector<double> cuts (getCuts ( parser.getOptArg("-cuts") ) );
	double xmin,xmax,dx;
	unsigned int i;

	// Use matlabformat?
	bool matlabformat ( parser.getOptSet ( "--matlab" ) );


	// Output stuff
	std::vector< std::vector<double> > predicted;
	std::vector<double> devianceresiduals;

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
		std::cerr << "\nparams:  ";
		for (i=0; i<theta.size(); i++) std::cerr << theta[i] << " ";
		std::cerr << "\n";
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

		// Replace cuts with the derived thresholds at the cuts
		for ( i=0; i<cuts.size(); i++ ) {
			cuts[i] = pmf->getThres ( theta, cuts[i] );
		}

		if ( verbose ) std::cerr << "\n";

		// Print output
		predicted = std::vector< std::vector<double> > ( data->getNblocks(), std::vector<double>(2) );
		xmin = 1e10; xmax=-1e10;
		for ( i=0; i<data->getNblocks(); i++ ) {
			predicted[i][0] = data->getIntensity(i);
			predicted[i][1] = pmf->evaluate ( data->getIntensity(i), theta );
			if ( data->getIntensity(i) < xmin ) xmin = data->getIntensity(i);
			if ( data->getIntensity(i) > xmax ) xmax = data->getIntensity(i);
		}
		print ( predicted, matlabformat, "prediction", ofile );
		predicted = std::vector< std::vector<double> > ( 100, std::vector<double>(2) );
		dx = xmax-xmin;
		dx /= 99;
		for ( i=0; i<100; i++ ) {
			predicted[i][0] = xmin + i*dx;
			predicted[i][1] = pmf->evaluate ( predicted[i][0], theta );
		}
		print ( predicted, matlabformat, "pmf", ofile );

		devianceresiduals = pmf->getDevianceResiduals ( theta, data );
		print ( devianceresiduals,                              matlabformat, "devianceresiduals", ofile );
		print ( cuts,                                           matlabformat, "thres", ofile );
		print ( pmf->deviance ( theta, data ),                  matlabformat, "deviance", ofile );
		print ( pmf->getRpd ( devianceresiduals, theta, data ), matlabformat, "rpd", ofile );
		print ( pmf->getRkd ( devianceresiduals, data ),        matlabformat, "rkd", ofile );

		fname = parser.popArg();

		// Clean up
		delete data;
		delete pmf;
	}

	return 0;
}
