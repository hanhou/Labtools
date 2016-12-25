/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#include <iostream>
#include <cstdlib>
#include "psychometric.h"
#include "mclist.h"
#include "bootstrap.h"
#include "testing.h"
#include "mcmc.h"
#include "getstart.h"
#include "integrate.h"

#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <signal.h>


int PsychometricValues ( TestSuite* T ) {
	int failures(0),i,j;
	char message[40];
	std::vector <double> x ( 6 );
	std::vector <int>    n ( 6, 50 );
	std::vector <int>    k ( 6 );
	std::vector <double> p ( 6 );

	std::vector <double> prm(3);
	std::vector <double> dl ( 3 );
	double d,l;
	std::vector<double> prm1 ( 4 );
	std::vector<double> dl1 ( 4 );

	// Set up data
	x[0] =  0.; x[1] =  2.; x[2] =  4.; x[3] =  6.; x[4] =  8.; x[5] = 10.;
	k[0] = 29;  k[1] = 31;  k[2] = 36;  k[3] = 42;  k[4] = 46;  k[5] = 49;
	p[0] = 0.5311852; p[1] = 0.60013209; p[2] = 0.74;
	p[3] = 0.87986791;p[4] = 0.9488148;  p[5] = 0.97136662;
	PsiData * data = new PsiData (x,n,k,2);

	// Set up psychometric function
	abCore * core = new abCore();
	PsiLogistic * sigmoid = new PsiLogistic();
	PsiPsychometric * pmf = new PsiPsychometric ( 2, core, sigmoid );
	prm[0] = 4; prm[1] = 1.5; prm[2] = 0.02;

	// Test forward probabilities
	for ( i=0; i<6; i++ ) {
		sprintf(message,"PsychometricValues at x=%g", x[i]);
		failures += T->isequal ( pmf->evaluate(x[i],prm), p[i], message);
	}

	// Test likelihood
	failures += T->isequal ( pmf->negllikeli(prm,data), 11.996474658154325, "PsychometricValues likelihood");
	
	// Test likelihood gradient
	dl = pmf->dnegllikeli ( prm, data );
	l  = pmf->negllikeli  ( prm, data );
	for ( i=0; i<3; i++ ) {
		prm[i] += 1e-5;
		d = pmf->negllikeli ( prm, data ) - l;
		d /= 1e-5;
		prm[i] -= 1e-5;
		failures += T->isequal ( dl[i], d, "PsychometricValues likelihood derivative", .05 );
	}

	// Test likelihood hessian
	Matrix *H = pmf->ddnegllikeli ( prm, data );
	dl = pmf->dnegllikeli ( prm, data );
	for ( i=0; i<3; i++ ) {
		prm[i] += 1e-9;
		dl1 = pmf->dnegllikeli ( prm, data );
		prm[i] -= 1e-9;
		for ( j=0; j<3; j++ ) {
			d = dl1[j] - dl[j];
			d /= 1e-9;
			failures += T->isequal ( (*H)(i,j), -d, "Psychometric function likelihood Hessian", .1 );
		}
	}
	delete H;
	delete pmf;

	// Yes no task
	pmf = new PsiPsychometric ( 1, core, sigmoid );
	prm1[0] = 4; prm1[1] = 1.5; prm1[2] = 0.02; prm1[3] = 0.5;

	// Test likelihood
	failures += T->isequal ( pmf->negllikeli(prm1,data), 11.996474658154325, "PsychometricValues likelihood-1afc");

	// Test likelihood gradient
	dl1 = pmf->dnegllikeli ( prm1, data );
	l  = pmf->negllikeli  ( prm1, data );
	for ( i=0; i<4; i++ ) {
		prm1[i] += 1e-5;
		d = pmf->negllikeli ( prm1, data ) - l;
		d /= 1e-5;
		prm1[i] -= 1e-5;
		failures += T->isequal ( dl1[i], d, "PsychometricValues likelihood-1afc derivative", .05 );
	}
	delete pmf;

	pmf = new BetaPsychometric ( 2, core, sigmoid );

	std::vector<double> bprm(4);
	bprm[0] = 4; bprm[1] = 1.5; bprm[2] = .02; bprm[3] = 1;

	// Observe, that beta likelihood can also be > 1 implying that both signs for log likelihood are possible
	failures += T->isequal ( pmf->negllikeli(bprm,data), -11.3918, "PsychometricValues beta likelihood", 1e-4);

	// Test likelihood gradient
	dl = pmf->dnegllikeli ( bprm, data );
	l  = pmf->negllikeli ( bprm, data );
	for ( i=0; i<4; i++ ) {
		bprm[i] += 1e-5;
		d = pmf->negllikeli ( bprm, data ) - l;
		d /= 1e-5;
		bprm[i] -= 1e-5;
		failures += T->isequal ( dl[i], d, "PsychometricValues beta likelihood derivative", .05 );
	}

	// test likelihood hessian
	H = pmf->ddnegllikeli ( bprm, data );
	// H->print();
	for ( i=0; i<4; i++ ) {
		bprm[i] += 1e-9;
		dl1 = pmf->dnegllikeli ( bprm, data );
		bprm[i] -= 1e-9;
		for ( j=0; j<4; j++ ) {
			d = dl1[j] - dl[j];
			d /= 1e-9;
			// failures += T->isequal ( log10((*H)(i,j)/ d), 0, "Psychometric Values beta likelihood Hessian", .12 );
			failures += T->isequal_rel ( (*H)(i,j), d, "Psychometric Values beta likelihood Hessian", .12 );
		}
	}
	delete H;

	failures += T->ismore ( pmf->deviance ( bprm, data ), 0, "Psychometric Values beta deviance" );
	delete pmf;

	delete core;
	delete sigmoid;
	delete data;

	return failures;
}

int DerivativeCheck ( TestSuite * T ) {
	/* Check all derivatives */
	int failures (0);
	unsigned int i,j, coreindex;
	double d;
	double y1,y0;
	double x(2.);
	std::vector<double> prm (3);
	prm[0] = 4.; prm[1] = 1.5; prm[2] = .02;
	char msg[500];
	std::vector <double> intensity ( 6 );
	std::vector <int>    n ( 6, 50 );
	std::vector <int>    k ( 6 );
	intensity[0] =  0.1; intensity[1] =  2.; intensity[2] =  4.; intensity[3] =  6.; intensity[4] =  8.; intensity[5] = 10.;
	k[0] = 29;  k[1] = 31;  k[2] = 36;  k[3] = 42;  k[4] = 46;  k[5] = 49;
	PsiData * data = new PsiData (intensity,n,k,2);

	// Cores
	PsiCore* core;
	std::vector<PsiCore*> cores ( 8 );
	std::vector<char*>    corenames ( 8 );
	cores[0] = new abCore (data);        corenames[0] = new char [20]; sprintf ( corenames[0], "abCore");
	cores[1] = new linearCore (data);    corenames[1] = new char [20]; sprintf ( corenames[1], "linearCore");
	cores[2] = new logCore (data);       corenames[2] = new char [20]; sprintf ( corenames[2], "logCore");
	cores[3] = new mwCore (data);        corenames[3] = new char [20]; sprintf ( corenames[3], "mwCore");
	cores[4] = new polyCore (data);      corenames[4] = new char [20]; sprintf ( corenames[4], "polyCore");
	cores[5] = new weibullCore (data);   corenames[5] = new char [20]; sprintf ( corenames[5], "weibullCore");
	cores[6] = new NakaRushton ( data ); corenames[6] = new char [20]; sprintf ( corenames[6], "NakaRushton");
	cores[7] = new mwCore (data,3);      corenames[3] = new char [20]; sprintf ( corenames[3], "mwCore (GumbelR)");

	for ( coreindex=0; coreindex<cores.size(); coreindex++ ) {
		core = cores[coreindex];
		// First derivative
		y0 = core->g ( x, prm );
		for ( i=0; i<3; i++ ) {
			prm[i] += 1e-7;
			y1 = core->g ( x, prm );
			prm[i] -= 1e-7;
			d = y1-y0; d /= 1e-7;
			sprintf ( msg, "%s 1st derivative w.r.t. prm %d", corenames[coreindex], i );
			failures += T->isequal ( core->dg ( x, prm, i ), d, msg, 1e-3 );
		}
		// Second derivative
		for ( i=0; i<3; i++ ) {
			y0 = core->dg ( x, prm, i );
			for ( j=0; j<3; j++ ) {
				prm[j] += 1e-7;
				y1 = core->dg ( x, prm, i );
				prm[j] -= 1e-7;
				d = y1-y0; d /= 1e-7;
				sprintf ( msg, "%s 2nd derivative w.r.t. prm %d and %d", corenames[coreindex], i, j );
				failures += T->isequal ( core->ddg ( x, prm, i, j ), d, msg, 2.5*1e-3 );
			}
		}
		delete cores[coreindex];
		delete corenames[coreindex];
	}

	// Sigmoids
	PsiSigmoid * sigmoid;
	std::vector<PsiSigmoid*> sigmoids ( 6 );
	std::vector<char*>       sigmoidnames ( 6 );
	sigmoids[0] = new PsiCauchy ();      sigmoidnames[0] = new char [20]; sprintf ( sigmoidnames[0], "PsiCauchy" );
	sigmoids[1] = new PsiExponential (); sigmoidnames[1] = new char [20]; sprintf ( sigmoidnames[1], "PsiExponential" );
	sigmoids[2] = new PsiGauss ();       sigmoidnames[2] = new char [20]; sprintf ( sigmoidnames[2], "PsiGauss" );
	sigmoids[3] = new PsiGumbelL ();     sigmoidnames[3] = new char [20]; sprintf ( sigmoidnames[3], "PsiGumbelL" );
	sigmoids[4] = new PsiGumbelR ();     sigmoidnames[4] = new char [20]; sprintf ( sigmoidnames[4], "PsiGumbelR" );
	sigmoids[5] = new PsiLogistic ();    sigmoidnames[5] = new char [20]; sprintf ( sigmoidnames[5], "PsiLogistic" );
	sigmoids[6] = new PsiId ();          sigmoidnames[6] = new char [20]; sprintf ( sigmoidnames[6], "PsiId" );

	for ( coreindex=0; coreindex<sigmoids.size(); coreindex++ ) {
		sigmoid = sigmoids[coreindex];
		// First derivative
		y0 = sigmoid->f ( x );
		y1 = sigmoid->f ( x+1e-7 );
		d = y1-y0; d /= 1e-7;
		sprintf ( msg, "%s 1st derivative", sigmoidnames[coreindex] );
		failures += T->isequal ( sigmoid->df ( x ), d, msg, 1e-3 );
		// Second derivative
		y0 = sigmoid->df ( x );
		y1 = sigmoid->df ( x+1e-7 );
		d = y1-y0; d /= 1e-7;
		sprintf ( msg, "%s 2nd derivative", sigmoidnames[coreindex] );
		failures += T->isequal ( sigmoid->ddf ( x ), d, msg, 1e-3 );
		delete sigmoids[coreindex];
		delete sigmoidnames[coreindex];
	}

	// Psychometric function
	core = new abCore(); sigmoid = new PsiLogistic ();
	PsiPsychometric * pmf = new PsiPsychometric ( 2, core, sigmoid );
	delete core;
	delete sigmoid;
	// First derivatives
	y1 = pmf->evaluate ( x, prm );
	for ( i=0; i<3; i++ ) {
		prm[i] += 1e-7;
		y0 = pmf->evaluate ( x, prm );
		prm[i] -= 1e-7;
		d = y0-y1; d /= 1e-7;
		sprintf ( msg, "Psychometric function 1st derivative w.r.t. prm %d", i );
		failures += T->isequal ( pmf->dpredict ( prm, x, i ), d, msg );
	}
	// Second derivatives
	for ( i=0; i<3; i++ ) {
		for ( j=0; j<3; j++ ) {
			y0 = pmf->dpredict ( prm, x, i );
			prm[j] += 1e-7;
			y1 = pmf->dpredict ( prm, x, i );
			prm[j] -= 1e-7;
			d = y1 - y0; d /= 1e-7;
			sprintf ( msg, "Psychometric function 2nd derivative w.r.t. prm %d and %d", i, j );
			failures += T->isequal ( pmf->ddpredict ( prm, x, i, j ), d, msg );
		}
	}
	delete pmf;

	// Special functions
	// psi = d log(Gamma)/ d x
	for ( x=.5; x<55; x+=5 ) {
		y1 = gammaln ( x );
		y0 = gammaln ( x+1e-7 );
		d = y0-y1; d /= 1e-7;
		sprintf ( msg, "psi function at x=%g", x );
		failures += T->isequal ( psi(x), d, msg, 1e-5 );
	}

	// digamma = d psi / dx
	for ( x=.5; x<55; x+=5 ) {
		y1 = psi ( x );
		y0 = psi ( x+1e-7 );
		d = y0-y1; d /= 1e-7;
		sprintf ( msg, "digamma function at x=%g", x );
		failures += T->isequal ( digamma(x), d, msg, 1e-5 );
	}

	return failures;
}

int BetaModelTest ( TestSuite * T ) {
	int failures ( 0 );

	std::vector<double> x ( 6 );
	std::vector<int>    n ( 6, 50 );
	std::vector<int>    k ( 6 );

	// Set up data
	x[0] =  0.; x[1] =  2.; x[2] =  4.; x[3] =  6.; x[4] =  8.; x[5] = 10.;
	k[0] = 24;  k[1] = 32;  k[2] = 40;  k[3] = 48;  k[4] = 50;  k[5] = 48;
	PsiData * data = new PsiData (x,n,k,2);

	PsiCore * core = new mwCore ();
	PsiSigmoid * sigmoid = new PsiLogistic ();

	PsiPrior * prior = new BetaPrior ( 2, 30 );
	PsiPsychometric * pmf = new BetaPsychometric ( 2, core, sigmoid );
	std::vector<double> prm(4);
	prm[0] = 4; prm[1] = 0.8; prm[2] = 0.02; prm[3] = .99;
	pmf->setPrior( 2, prior );
	delete prior;
	prior = new UniformPrior ( 0, 1 );
	pmf->setPrior( 3, prior );
	delete prior;

	PsiOptimizer optimizer ( pmf, data );
	std::vector<double> mapestimate = optimizer.optimize ( pmf, data );

	failures += T->isequal ( mapestimate[0], 3.99318,   "Beta model MAP estimate m"     , 1e-5 );
	failures += T->isequal ( mapestimate[1], 3.87268,   "Beta model MAP estimate w"     , 1e-5 );
	failures += T->isequal ( mapestimate[2], 0.02,      "Beta model MAP estimate lambda", 1e-5 );
	failures += T->isequal ( mapestimate[3], 0.99999,  "Beta model MAP estimate nu"    , 1e-5 );

	setSeed(0);
	GenericMetropolis * gmS = new GenericMetropolis ( pmf, data, new GaussRandom() );
	gmS->setTheta ( prm );

	MCMCList pilot ( gmS->sample(1000) );
	gmS->findOptimalStepwidth(pilot);
	MCMCList post = gmS->sample(3000);

	failures += T->isequal ( post.getMean ( 0 ), 3.33222,   "Beta model MEAN estimate m"     , 1e-5 );
	failures += T->isequal ( post.getMean ( 1 ), 4.31489,   "Beta model MEAN estimate w"     , 1e-5 );
	failures += T->isequal ( post.getMean ( 2 ), 0.0395864, "Beta model MEAN estimate lambda", 1e-5 );
	failures += T->isequal ( post.getMean ( 3 ), 0.622763,  "Beta model MEAN estimate nu"    , 1e-5 );

	return failures;
}

int OptimizerSolution ( TestSuite * T ) {
	int failures(0);
	unsigned int i;
	double deviance(0);

	/************************
	 * 2AFC
	 */
	std::vector <double> x ( 6 );
	std::vector <int>    n ( 6, 50 );
	std::vector <int>    k ( 6 );

	// Set up data
	x[0] =  0.; x[1] =  2.; x[2] =  4.; x[3] =  6.; x[4] =  8.; x[5] = 10.;
	k[0] = 24;  k[1] = 32;  k[2] = 40;  k[3] = 48;  k[4] = 50;  k[5] = 48;
	PsiData * data = new PsiData (x,n,k,2);

	// Set up psychometric function
	abCore * core = new abCore();
	PsiLogistic * sigmoid = new PsiLogistic();
	PsiPsychometric * pmf = new PsiPsychometric ( 2, core, sigmoid );
	std::vector<double> prm(4);
	prm[0] = 4; prm[1] = 0.8; prm[2] = 0.02;
	PsiPrior *prior = new UniformPrior(0.,0.1);
	pmf->setPrior( 2, prior );

	// Optimizer
	PsiOptimizer *opt = new PsiOptimizer (pmf, data);

	std::vector<double> solution (4);
	solution = opt->optimize(pmf,data);
	std::vector<double> devianceresiduals (pmf->getDevianceResiduals ( solution, data ));
	for ( i=0; i<devianceresiduals.size(); i++ ) {
		deviance += devianceresiduals[i]*devianceresiduals[i];
	}

	// Test against old psignifit:
	// tests/testdata_from_psignifit$ psignifit testprefs test.dat
	failures += T->isequal(solution[0],3.29515,"OptimizerSolution 2AFC alpha",1e-2);
	failures += T->isequal(solution[1],0.960018,"OptimizerSolution 2AFC beta",1e-2);
	failures += T->isequal(solution[2],0.0190379,"OptimizerSolution 2AFC lambda",1e-3);

	failures += T->isequal(pmf->deviance(solution,data),3.98476,"OptimizerSolution 2AFC deviance",1e-2);
	failures += T->isequal(pmf->deviance(solution,data),deviance,"OptimizerSolution 2AFC deviance sum", 1e-7);

	failures += T->isequal(pmf->getRpd(devianceresiduals,solution,data),0.155395,"OptimizerSolution 2AFC Rpd",1e-2);
	failures += T->isequal(pmf->getRkd(devianceresiduals,data),-0.320889,"OptimizerSolution 2AFC Rkd",1e-2);

	delete pmf;
	delete data;
	delete opt;

	/************************
	 * Yes/No
	 */
	std::clog << "\n";

	// Set up data
	x[0] =  0.; x[1] =  2.; x[2] =  4.; x[3] =  6.; x[4] =  8.; x[5] = 10.;
	k[0] = 3;  k[1] = 10;  k[2] = 34;  k[3] = 45;  k[4] = 50;  k[5] = 50;
	data = new PsiData (x,n,k,2);

	// Set up psychometric function
	pmf = new PsiPsychometric ( 1, core, sigmoid );
	prm[0] = 4; prm[1] = 0.8; prm[2] = 0.02; prm[3] = 0.1;
	pmf->setPrior( 2, prior);
	pmf->setPrior( 3, prior);

	// Optimizer
	opt = new PsiOptimizer (pmf, data);

	solution = opt->optimize(pmf,data);

	devianceresiduals = pmf->getDevianceResiduals ( solution, data );
	deviance = 0;
	for ( i=0; i<devianceresiduals.size(); i++ ) {
		deviance += devianceresiduals[i]*devianceresiduals[i];
	}

	// Check against classical psignifit:
	// tests/testdata_from_psignifit$ psignifit testprefs_yn test_yn.dat
	failures += T->isequal(solution[0],3.43942,"OptimizerSolution Y/N alpha",1e-2);
	failures += T->isequal(solution[1],0.988357,"OptimizerSolution Y/N beta", 5*1e-3);
	failures += T->isequal(solution[2],3.61604e-8,"OptimizerSolution Y/N lambda",5*1e-3);
	failures += T->isequal(solution[3],0.028935,"OptimizerSolution Y/N gamma",5*1e-3);

	failures += T->isequal(pmf->deviance(solution,data),2.08172,"OptimizerSolution Y/N deviance",1e-2);
	failures += T->isequal(pmf->deviance(solution,data),deviance,"OptimizerSolution Y/N deviance sum", 1e-7);

	failures += T->isequal(pmf->getRpd(devianceresiduals,solution,data),0.217146,"OptimizerSolution Y/N Rpd",1e-2);
	failures += T->isequal(pmf->getRkd(devianceresiduals,data),-0.477967,"OptimizerSolution Y/N Rkd",2e-2);

	delete pmf;
	delete opt;

	// Yes/No with gamma==lambda
	std::cerr << "\n";
	pmf = new PsiPsychometric ( 1, core, sigmoid );
	pmf->setgammatolambda ();
	opt = new PsiOptimizer ( pmf, data );
	pmf->setPrior( 2, prior);
	solution = opt->optimize(pmf,data);
	
	failures += T->isequal ( solution[0], 3.3044829936747959, "Optimizer Solution Yes/No gamma=lambda, alpha", 1e-3 );
	failures += T->isequal ( solution[1], 1.0698582800005716, "Optimizer Solution Yes/No gamma=lambda, beta", 1e-3 );
	failures += T->isequal ( solution[2], 4.4302035660369048e-10, "Optimizer Solution Yes/No gamma=lambda, lambda", 1e-12 );
	failures += T->isequal (pmf->deviance(solution,data),2.49387706984, "Optimizer Solution Yes/No gamma=lambda, deviance", 1e-2);

	delete pmf;
	delete data;
	delete opt;
	delete core;
	delete sigmoid;
	delete prior;

	return failures;
}

int InitialParametersTest ( TestSuite * T ) {
	int i, failures(0);
	std::vector<double> x ( 4 );
	std::vector<int>    k ( 4 );
	std::vector<int>    n ( 4, 20 );
	PsiPsychometric * pmf;
	PsiData * data;
	std::vector<double> prm ( 3 );

	for (i=1; i<=4; i++) {
		x[i-1] = i;
		k[i-1] = 10+(i-1)*3;
	}
	data = new PsiData ( x, n, k, 2 );

	abCore * core = new abCore();
	PsiLogistic * sigmoid = new PsiLogistic();
	pmf = new PsiPsychometric ( 2, core, sigmoid );
	prm = pmf->getStart ( data );
	failures += T->ismore ( prm[1], 0, "PsiPsychometric->getStart() for increasing data" );

	delete data;

	for (i=1; i<=4; i++) {
		x[i-1] = 5-i;
	}
	data = new PsiData ( x, n, k, 2 );

	prm = pmf->getStart ( data );
	failures += T->isless ( prm[1], 0, "PsiPsychometric->getStart() for decreasing data" );

	delete data;
	delete pmf;
	delete core;
	delete sigmoid;

	return failures;
}

int BootstrapTest ( TestSuite * T ) {
	setSeed ( 0 );
	int failures(0);
	unsigned int i;
	std::vector<double> x ( 6 );
	std::vector<int>    n ( 6, 50 );
	std::vector<int>    k ( 6 );
	char testname[40];

	// Set up data
	x[0] =  0.; x[1] =  2.; x[2] =  4.; x[3] =  6.; x[4] =  8.; x[5] = 10.;
	k[0] = 24;  k[1] = 32;  k[2] = 40;  k[3] = 48;  k[4] = 50;  k[5] = 48;
	PsiData * data = new PsiData (x,n,k,2);

	// Set up psychometric function
	abCore * core = new abCore();
	PsiLogistic * sigmoid = new PsiLogistic();
	PsiPrior * prior = new UniformPrior ( 0, .1 );
	PsiPsychometric * pmf = new PsiPsychometric ( 2, core, sigmoid );
	std::vector<double> prm(4);
	prm[0] = 4; prm[1] = 0.8; prm[2] = 0.02;
	pmf->setPrior( 2, prior );

	std::vector<double> cuts (1, 0.5);
	// BootstrapList boots = bootstrap ( 9999, data, pmf, cuts );
	BootstrapList boots = bootstrap ( 999, data, pmf, cuts );

	// Check against psignifit results
	// These values are subject to statistical variation. "equality" is defined relatively coarse
	failures += T->isless(boots.getAcc_t(0),     0.018662,"Acceleration constant (threshold)");
	failures += T->isequal(boots.getBias_t(0), -0.0928786,"Bias (threshold)",            .01);
	failures += T->isequal(boots.getThres(.1,0), 2.65266,"th(.1)",                        .05);
	failures += T->isequal(boots.getThres(.9,0), 3.89757,"th(.9)",                        .05);

	failures += T->isequal(boots.getAcc_s(0),    -0.000155314, "Acceleration constant (slope)", .01);
	failures += T->isequal(boots.getBias_s(0),    -0.0501536,   "Bias (slope)",                  .01);
	failures += T->isequal(boots.getSlope(0.1,0), 0.181289,    "sl(.1)",                        .01);
	failures += T->isequal(boots.getSlope(0.9,0), 0.497512,    "sl(.9)",                        .01);

	failures += T->isequal(boots.getDeviancePercentile(0.975),9.67016,"Deviance limits",.5);
	failures += T->isequal(boots.percRpd(.025), -0.451653, "Rpd( 2.5%)", .1); // Testing mean and standard error
	failures += T->isequal(boots.percRpd(.975), 0.632072, "Rpd(97.5%)",  .1);
	failures += T->isequal(boots.percRkd(.025), -0.932597, "Rkd( 2.5%)", .1);
	failures += T->isequal(boots.percRkd(.975), 0.601175, "Rkd(97.5%)",  .1);

	// Check for influential observations and outliers
	JackKnifeList jackknife = jackknifedata (data, pmf);

	std::vector<double> ci_lower ( pmf->getNparams() ), ci_upper ( pmf->getNparams() );
	for ( i=0; i<pmf->getNparams(); i++ ) {
		ci_lower[i] = boots.getPercentile(0.025,i);
		ci_upper[i] = boots.getPercentile(0.975,i);
	}

	for ( i=0; i<6; i++ ) {
		sprintf(testname,"influential %d",i);
		failures += T->conditional(jackknife.influential(i,ci_lower,ci_upper)<1,testname);
		sprintf(testname,"outliers %d",i);
		failures += T->conditional(!jackknife.outlier(i),testname);
	}

	delete core;
	delete sigmoid;
	delete prior;
	delete pmf;
	delete data;

	return failures;
}

int MCMCTest ( TestSuite * T ) {
	int failures ( 0 );

	std::vector<double> x ( 6 );
	std::vector<int>    n ( 6, 50 );
	std::vector<int>    k ( 6 );

	// Set up data
	x[0] =  0.; x[1] =  2.; x[2] =  4.; x[3] =  6.; x[4] =  8.; x[5] = 10.;
	k[0] = 24;  k[1] = 32;  k[2] = 40;  k[3] = 48;  k[4] = 50;  k[5] = 48;
	PsiData * data = new PsiData (x,n,k,2);

	// Set up psychometric function
	PsiCore * core = new abCore ();
	PsiSigmoid * sigmoid = new PsiLogistic();
	PsiPrior * prior = new UniformPrior ( 0, .1 );
	PsiPsychometric * pmf = new PsiPsychometric ( 2, core, sigmoid );
	std::vector<double> prm(3);
	prm[0] = 4; prm[1] = 0.8; prm[2] = 0.02;
	pmf->setPrior( 2, prior );

	MetropolisHastings * mhS = new MetropolisHastings( pmf, data, new GaussRandom() );
	mhS->setTheta( prm );
	mhS->setStepSize(0.1,0);
	mhS->setStepSize(0.1,1);
	mhS->setStepSize(0.001,2);

	GenericMetropolis * gmS = new GenericMetropolis ( pmf, data, new GaussRandom() );
	gmS->setTheta ( prm );

	HybridMCMC * S = new HybridMCMC ( pmf, data, 20 );
	S->setTheta ( prm );
	// This gives rather bad sampling but at least it gives something
	// We don't use the HybridMCMC anyhow
	S->setStepSize ( 0.013, 0 );
	S->setStepSize ( 0.007, 1 );
	S->setStepSize ( 0.001, 2 );

	srand48(0);
	MCMCList post ( S->sample(1000) );
	srand48(0);
	MCMCList mhpost ( mhS->sample(1000) );
	srand48(0);
	MCMCList pilot ( mhS->sample(1000) );
	gmS->findOptimalStepwidth(pilot);
	srand48(0);
	MCMCList gmpost = gmS->sample(1000);

	/*
	int i,j;
	for ( i=0; i<1000; i++ ) {
		for ( j=0; j<3; j++)
			std::cout << " " << post.getEst ( i, j );
		std::cout << "\n";
	}
	*/

	failures += T->isequal ( post.getMean(0), 3.58027, "Hybrid MCMC alpha", .3 );
	failures += T->isequal ( post.getMean(1), 0.909616, "Hybrid MCMC beta", .2 );
	failures += T->isequal ( post.getMean(2), 0.0217217, "Hybrid MCMC lambda", .02 );
	failures += T->isequal ( mhpost.getMean(0), 3.22372, "Metropolis Hastings alpha", .2 );
	failures += T->isequal ( mhpost.getMean(1), 1.12734, "Metropolis Hastings beta", .2 );
	failures += T->isequal ( mhpost.getMean(2), 0.0199668, "Metropolis Hastings lambda", .02 );
	failures += T->isequal ( gmpost.getMean(0), 3.22372, "Generic Metropolis MCMC alpha", .2 );
	failures += T->isequal ( gmpost.getMean(1), 1.12734, "Generic Metropolis MCMC beta", .2 );
	failures += T->isequal ( gmpost.getMean(2), 0.0199668, "Generic Metropolis MCMC lambda", .02 );

	delete core;
	delete sigmoid;
	delete prior;
	delete pmf;
	delete mhS;
	delete S;
	delete data;

	return failures;
}

int PriorTest ( TestSuite * T ) {
	int failures ( 0 );
	PsiPrior * prior;

	prior = new PsiPrior;
	failures += T->isequal ( prior->pdf ( 0 ), 1, "Flat prior at 0" );
	failures += T->isequal ( prior->dpdf ( 0 ) , 0, "Flat prior derivative at 0" );
	delete prior;

	prior = new UniformPrior ( 0, 1 );
	failures += T->isequal ( prior->pdf ( -.5 ) , 0, "Uniform prior at -0.5" );
	failures += T->isequal ( prior->pdf ( .5 ) , 1,  "Uniform prior at 0.5" );
	failures += T->isequal ( prior->pdf ( 1.5 ) , 0, "Uniform prior at 1.5" );
	failures += T->isequal ( prior->dpdf ( -.5 ) , 0, "Uniform prior derivative at -0.5" );
	failures += T->isequal ( prior->dpdf ( .5 ) , 0,  "Uniform prior derivative at 0.5" );
	failures += T->isequal ( prior->dpdf ( 1.5 ) , 0, "Uniform prior derivative at 1.5" );
	delete prior;

	prior = new GaussPrior ( 0, 1 );
	failures += T->isequal ( prior->pdf ( -1 ), 0.24197072, "Gaussian prior at -1" );
	failures += T->isequal ( prior->pdf ( 0 ), 0.39894228, "Gaussian prior at 0" );
	failures += T->isequal ( prior->pdf ( 1 ), 0.24197072, "Gaussian prior at 1" );
	failures += T->isequal ( prior->dpdf ( -1 ), 0.24197072, "Gaussian prior derivative at -1" );
	failures += T->isequal ( prior->dpdf ( 0 ), 0, "Gaussian prior derivative at 0" );
	failures += T->isequal ( prior->dpdf ( 1 ), -0.24197072, "Gaussian prior derivative at 1" );
	delete prior;

	prior = new BetaPrior ( 1.5, 3. );
	failures += T->isequal ( prior->pdf ( -.1 ), 0, "BetaPrior at 0" );
	failures += T->isequal ( prior->pdf ( .1 ), 1.68094822, "BetaPrior at 0.1" );
	failures += T->isequal ( prior->pdf ( .5 ), 1.16009706, "BetaPrior at 0.5" );
	failures += T->isequal ( prior->pdf ( 1.1 ), 0, "BetaPrior at 1.1" );
	failures += T->isequal ( prior->dpdf ( -.1 ), 0, "BetaPrior derivative at 0" );
	failures += T->isequal ( prior->dpdf ( .1 ), 12.14018158, "BetaPrior derivative at 0.1" );
	failures += T->isequal ( prior->dpdf ( .5 ), 5.80048531, "BetaPrior derivative at 0.5" );
	failures += T->isequal ( prior->dpdf ( 1.1 ), 0, "BetaPrior derivative at 1.1" );
	delete prior;

	prior = new GammaPrior ( 1.5, 3. );
	failures += T->isequal ( prior->pdf ( -0.5 ), 0., "GammaPrior at -0.5" );
	failures += T->isequal ( prior->pdf ( 0.5 ), 0.12997977, "GammaPrior at 0.5" );
	failures += T->isequal ( prior->pdf ( 1.0 ), 0.15559955, "GammaPrior at 1.0" );
	failures += T->isequal ( prior->pdf ( 1.5 ), 0.16131382, "GammaPrior at 1.5" );
	failures += T->isequal ( prior->dpdf ( -0.5 ), 0., "GammaPrior derivative at -0.5" );
	failures += T->isequal ( prior->dpdf ( 0.5 ), 0.08665318, "GammaPrior derivative at 0.5" );
	failures += T->isequal ( prior->dpdf ( 1.0 ), 0.02593326, "GammaPrior derivative at 1.0" );
	failures += T->isequal ( prior->dpdf ( 1.5 ), 0., "GammaPrior derivative at 1.5" );
	delete prior;

	prior = new nGammaPrior ( 1.5, 3. );
	failures += T->isequal ( prior->pdf ( 0.5 ), 0., "nGammaPrior at 0.5" );
	failures += T->isequal ( prior->pdf ( -0.5 ), 0.12997977, "nGammaPrior at -0.5" );
	failures += T->isequal ( prior->pdf ( -1.0 ), 0.15559955, "nGammaPrior at -1.0" );
	failures += T->isequal ( prior->pdf ( -1.5 ), 0.16131382, "nGammaPrior at -1.5" );
	failures += T->isequal ( prior->dpdf ( 0.5 ), 0., "nGammaPrior derivative at 0.5" );
	failures += T->isequal ( prior->dpdf ( -0.5 ), -0.08665318, "nGammaPrior derivative at -0.5" );
	failures += T->isequal ( prior->dpdf ( -1.0 ), -0.02593326, "nGammaPrior derivative at -1.0" );
	failures += T->isequal ( prior->dpdf ( -1.5 ), 0., "nGammaPrior derivative at -1.5" );
	delete prior;


	return failures;
}

int SigmoidTests ( TestSuite * T ) {
	int failures(0);
	PsiSigmoid * sigmoid;

	// Check gaussian cdf themselves
	failures += T->isequal(Phi(0),.5,"Phi(0)",1e-5);
	failures += T->isequal(invPhi(0.5),0.,"invPhi(0.5)",1e-5);
	failures += T->isequal(invPhi(Phi(.3)),.3,"invPhi(Phi(0.3))",1e-5);
	failures += T->isequal(Phi(invPhi(0.3)),.3,"Phi(invPhi(0.3))",1e-5);

	sigmoid = new PsiLogistic ();
	// Check specific function values
	// f should be 0.5 at 0 and close to 0 resp. 1 at low resp. high values
	failures += T->isequal ( sigmoid->f ( 0), 0.5,  "PsiLogistic->f( 0)" );
	failures += T->isless  ( sigmoid->f (-3), 0.05, "PsiLogistic->f(-3)" );
	failures += T->ismore  ( sigmoid->f ( 3), 0.95, "PsiLogistic->f( 3)" );
	// Check symmetry
	failures += T->isequal ( sigmoid->f ( 3), 1-sigmoid->f(-3), "PsiLogistic( 3)-PsiLogistic(-3)" );
	// Check saturation
	// f should have low derivative at values far from 0
	failures += T->isless  ( sigmoid->df( 3), 0.05, "PsiLogistic->df(3)");
	failures += T->isless  ( sigmoid->df(-3), 0.05, "PsiLogistic->df(-3)");
	// Check monotonicity
	double mindf(1e20),x,df;
	for ( x=-5; x<5; x+=0.1 )
		if ( (df=sigmoid->df(x))<mindf )
			mindf = df;
	failures += T->ismore ( mindf, 0, "PsiLogistic monotonically increasing" );
		
	// Check convexity
	// if x>0, ddf<0
	// if x<0, ddf>0
	failures += T->isless  ( sigmoid->ddf(3), 0,    "PsiLogistic->ddf(3)");
	failures += T->ismore  ( sigmoid->ddf(-3),0,    "PsiLogistic->ddf(-3)");
	delete sigmoid;

	sigmoid = new PsiGauss ();
	// Check specific function values
	failures += T->isequal ( sigmoid->f ( 0), 0.5,  "PsiGauss->f( 0)" );
	failures += T->isless  ( sigmoid->f (-3), 0.01, "PsiGauss->f(-3)" );
	failures += T->ismore  ( sigmoid->f ( 3), 0.99, "PsiGauss->f( 3)" );
	// Check symmetry
	failures += T->isequal ( sigmoid->f ( 3), 1-sigmoid->f(-3), "PsiGauss( 3)-PsiGauss(-3)" );
	// Check monotonicity
	mindf = 1e20;
	for ( x=-5; x<5; x+=0.1 )
		if ( (df=sigmoid->df(x))<mindf )
			mindf = df;
	failures += T->ismore ( mindf, 0, "PsiGaussian monotonically increasing" );
	// Check saturation
	failures += T->isless  ( sigmoid->df( 3), 0.01, "PsiGauss->df(3)");
	failures += T->isless  ( sigmoid->df(-3), 0.01, "PsiGauss->df(-3)");
	// Check convexity
	failures += T->isless  ( sigmoid->ddf(3), 0,    "PsiGauss->ddf(3)");
	failures += T->ismore  ( sigmoid->ddf(-3),0,    "PsiGauss->ddf(-3)");
	delete sigmoid;

	sigmoid = new PsiGumbelL ();
	// Check specific function values
	failures += T->isequal ( sigmoid->f (0), 0.63212055882855767, "PsiGumbelL->f(0)");
	failures += T->isequal ( sigmoid->f (-3), .048568007099546562, "PsiGumbelL->f(0)");
	failures += T->isequal ( sigmoid->f (3), .99999999810782125, "PsiGumbelL->f(0)");
	// Check asymmetry
	failures += T->ismore ( sigmoid->f ( 3 ), 1-sigmoid->f( -3 ), "PsiGumbelL( 3 )-PsiGumbelL (-3 )" );
	// Check monotonicity
	mindf = 1e20;
	for (x=-5; x<5; x+=.1 )
		if ( (df=sigmoid->df(x))<mindf )
			mindf = df;
	failures += T->ismore ( mindf, 0, "PsiGumbelL monotonically increasing" );
	// Check saturation
	failures += T->isless  ( sigmoid->df( 3), 0.01, "PsiGumbelL->df(3)");
	failures += T->isless  ( sigmoid->df(-3), 0.05, "PsiGumbelL->df(-3)");
	// Check convexity
	failures += T->isless  ( sigmoid->ddf(3), 0,    "PsiGumbelL->ddf(3)");
	failures += T->ismore  ( sigmoid->ddf(-3),0,    "PsiGumbelL->ddf(-3)");
	delete sigmoid;

	sigmoid = new PsiGumbelR ();
	// Check specific function values
	failures += T->isequal ( sigmoid->f (0), 0.36787944117144233, "PsiGumbelR->f(0)");
	failures += T->isequal ( sigmoid->f (-3), 1.8921786948382924e-09, "PsiGumbelR->f(0)");
	failures += T->isequal ( sigmoid->f (3), .95143199290045344, "PsiGumbelR->f(0)");
	// Check asymmetry
	failures += T->ismore ( 1-sigmoid->f ( -3 ), sigmoid->f( 3 ), "PsiGumbelR( -3 )-PsiGumbelR (3 )" );
	// Check monotonicity
	mindf = 1e20;
	for (x=-5; x<5; x+=.1 )
		if ( (df=sigmoid->df(x))<mindf )
			mindf = df;
	failures += T->ismore ( mindf, 0, "PsiGumbelR monotonically increasing" );
	// Check saturation
	failures += T->isless  ( sigmoid->df( 3), 0.05, "PsiGumbelR->df(3)");
	failures += T->isless  ( sigmoid->df(-3), 0.01, "PsiGumbelR->df(-3)");
	// Check convexity
	failures += T->isless  ( sigmoid->ddf(3), 0,    "PsiGumbelR->ddf(3)");
	failures += T->ismore  ( sigmoid->ddf(-3),0,    "PsiGumbelR->ddf(-3)");
	delete sigmoid;

	sigmoid = new PsiCauchy ();
	// Check specific function values
	failures += T->isequal ( sigmoid->f ( 0 ), 0.5, "PsiCauchy->f(0)" );
	double zalpha = -2*tan(M_PI*(0.1-0.5));
	failures += T->isequal ( sigmoid->f(-0.5*zalpha), 0.1, "PsiCauchy->f(-z(0.1)*(-.5) )" );
	failures += T->isequal ( sigmoid->f(0.5*zalpha),  0.9, "PsiCauchy->f( z(0.9)*.5 )" );
	// Check monotonicity
	mindf = 1e20;
	for (x=-5; x<5; x+=.1 ) {
		if ( (df=sigmoid->df(x))<mindf )
			mindf = df;
	}
	failures += T->ismore ( mindf, 0, "PsiCauchy monotonically increasing" );
	// Check saturation
	failures += T->isless  ( sigmoid->df( 3), 0.05, "PsiCauchy->df(3)");
	failures += T->isless  ( sigmoid->df(-3), 0.05, "PsiCauchy->df(-3)");
	// Check convexity
	failures += T->isless  ( sigmoid->ddf(3), 0,    "PsiCauchy->ddf(3)");
	failures += T->ismore  ( sigmoid->ddf(-3),0,    "PsiCauchy->ddf(-3)");
	delete sigmoid;

	sigmoid = new PsiId ();
	// A number of values
	failures += T->isequal ( sigmoid->f ( -1. ), -1., "PsiId->f(-1)" );
	failures += T->isequal ( sigmoid->f (  0. ),  0., "PsiId->f(0)"  );
	failures += T->isequal ( sigmoid->f (  1. ),  1., "PsiId->f( 1)" );

	// df
	failures += T->isequal ( sigmoid->df ( -1. ), 1., "PsiId->df(-1)" );
	failures += T->isequal ( sigmoid->df (  0. ), 1., "PsiId->df( 0)" );
	failures += T->isequal ( sigmoid->df (  1. ), 1., "PsiId->df( 1)" );

	// ddf
	failures += T->isequal ( sigmoid->ddf ( -1. ), 0., "PsiId->ddf(-1)" );
	failures += T->isequal ( sigmoid->ddf (  0. ), 0., "PsiId->ddf( 0)" );
	failures += T->isequal ( sigmoid->ddf (  1. ), 0., "PsiId->ddf( 1)" );

	// inv
	failures += T->isequal ( sigmoid->inv ( sigmoid->f ( -1. ) ), -1., "PsiId->inv(PsiId->f(-1))" );
	failures += T->isequal ( sigmoid->inv ( sigmoid->f (  0. ) ),  0., "PsiId->inv(PsiId->f( 0))" );
	failures += T->isequal ( sigmoid->inv ( sigmoid->f (  1. ) ),  1., "PsiId->inv(PsiId->f( 1))" );

	return failures;
}

int CoreTests ( TestSuite * T ) {
	int failures(0);
	PsiCore * core;
	PsiData * data;
	std::vector<double> *x;
	std::vector<int> *k,*n;
	unsigned int i;
	double th;
	std::vector<double> prm(4,0);
	std::vector<double> prm2(4,0);

	core = new abCore;
	prm[0] = 3.;
	prm[1] = 2.;
	failures += T->isequal ( core->g(3.,prm),0,                        "abCore at threshold");
	failures += T->isequal ( core->dgx(3.,prm),1./prm[1],              "abCore derivative stimulus");
	failures += T->isequal ( core->dg(3.,prm,0),-1./prm[1],            "abCore derivative 0");
	failures += T->isequal ( core->dg(3.,prm,1),0,                     "abCore derivative 1");
	failures += T->isequal ( core->ddg(3.,prm,0,0),0,                  "abCore 2nd derivative 0,0");
	failures += T->isequal ( core->ddg(3.,prm,0,1),1./(prm[1]*prm[1]), "abCore 2nd derivative 0,1");
	failures += T->isequal ( core->ddg(3.,prm,1,1),0,                  "abCore 2nd derivative 1,1");
	failures += T->isequal ( core->g(core->inv(2.,prm),prm),2,         "abCore inversion g(inv(2))");
	failures += T->isequal ( core->inv(core->g(2.,prm),prm),2,         "abCore inversion inv(g(2))");
	failures += T->isequal ( core->dinv(2.,prm,0),1.,                  "abCore inversion dinv(2,0)");
	failures += T->isequal ( core->dinv(2.,prm,1),2.,                  "abCore inversion dinv(2,1)");
	// TODO: Transform tests
	delete core;

	core = new mwCore (NULL, 1,0.1);
	prm[0] = 3.;
	prm[1] = 2.;
	failures += T->isequal ( core->g(3.,prm),0,                  "mwCore at threshold");
	failures += T->isequal ( core->dgx(3.,prm),log(9.),   "mwCore derivative stimulus");
	failures += T->isequal ( core->dg(3.,prm,0),-log(9.),        "mwCore derivative 0");
	failures += T->isequal ( core->dg(3.,prm,1),0,               "mwCore derivative 1");
	failures += T->isequal ( core->ddg(3.,prm,0,0),0,            "mwCore 2nd derivative 0,0");
	failures += T->isequal ( core->ddg(3.,prm,0,1),0.5*log(9.),  "mwCore 2nd derivative 0,1");
	failures += T->isequal ( core->ddg(3.,prm,1,1),0,            "mwCore 2nd derivative 1,1");
	failures += T->isequal ( core->g(core->inv(2.,prm),prm),2,   "mwCore inversion g(inv(2))");
	failures += T->isequal ( core->inv(core->g(2.,prm),prm),2,   "mwCore inversion inv(g(2))");
	failures += T->isequal ( core->dinv(2.,prm,0),1.,            "mwCore inversion dinv(2,0)");
	failures += T->isequal ( core->dinv(2.,prm,1),1./log(9.),    "mwCore inversion dinv(2,1)");
	// TODO: Transform Tests

	delete core;

	// Tests for m and w with all sigmoids
	std::vector<PsiSigmoid*> sigmoids ( 6 );
	std::vector<char*>       sigmnames ( 6 );
	sigmoids[0] = new PsiLogistic();    sigmnames[0] = new char [20]; sprintf ( sigmnames[0], "Logistic" );
	sigmoids[1] = new PsiGauss();       sigmnames[1] = new char [20]; sprintf ( sigmnames[1], "Gauss" );
	sigmoids[2] = new PsiGumbelL();     sigmnames[2] = new char [20]; sprintf ( sigmnames[2], "GumbelL" );
	sigmoids[3] = new PsiCauchy();      sigmnames[3] = new char [20]; sprintf ( sigmnames[3], "Cauchy" );
	sigmoids[4] = new PsiExponential(); sigmnames[4] = new char [20]; sprintf ( sigmnames[4], "Exponential" );
	sigmoids[5] = new PsiGumbelR();     sigmnames[5] = new char [20]; sprintf ( sigmnames[5], "GumbelR" );
	prm[0] = 4;
	prm[1] = 2;
	char message[40];

	prm[0] = 3.;
	for (i=0; i<sigmoids.size(); i++) {
		core = new mwCore ( NULL, sigmoids[i]->getcode(), 0.1 );
		sprintf(message,"mwCore (m) for Psi%s", sigmnames[i]);
		failures += T->isequal ( sigmoids[i]->f(core->g(prm[0],prm)), 0.5, message );
		sprintf(message,"mwCore (w) for Psi%s", sigmnames[i]);
		failures += T->isequal ( core->inv(sigmoids[i]->inv(0.9),prm) - core->inv(sigmoids[i]->inv(0.1),prm), prm[1], message );
		delete core;
		delete sigmoids[i];
		delete sigmnames[i];
	}

	core = new linearCore;
	prm2 = prm;
	th = -2./3;
	failures += T->isequal ( core->g(th,prm), 0,                   "linearCore at threshold");
	failures += T->isequal ( core->dgx(3.,prm),prm[0],             "linearCore derivative stimulus");
	failures += T->isequal ( core->dg(th,prm,0), th,               "linearCore derivative(0) at threshold");
	failures += T->isequal ( core->dg(th,prm,1), 1.,               "linearCore derivative(1) at threshold");
	failures += T->isequal ( core->ddg(th,prm,0,0), 0,             "linearCore derivative(0,0) at threshold");
	failures += T->isequal ( core->ddg(th,prm,1,0), 0,             "linearCore derivative(0,0) at threshold");
	failures += T->isequal ( core->ddg(th,prm,1,1), 0,             "linearCore derivative(0,0) at threshold");
	failures += T->isequal ( core->ddg(th,prm,0,1), core->ddg(th,prm,1,0), "linearCore 2nd derivative symmetry at threshold");
	failures += T->isequal ( core->inv(0,prm),th,                  "linearCore inverse threshold");
	failures += T->isequal ( core->g(core->inv(2.,prm),prm), 2.,   "linearCore inverse g(inv(2))");
	failures += T->isequal ( core->inv(core->g(2.,prm),prm), 2.,   "linearCore inverse inv(g(2))");
	prm2[0] += 1e-8;
	failures += T->isequal ( core->dinv(2.,prm,0), (core->inv(2,prm2)-core->inv(2.,prm))/1e-8,              "linearCore inverse derivative(0)");
	prm2[0] = prm[0];
	prm2[1] += 1e-8;
	failures += T->isequal ( core->dinv(2.,prm,1), (core->inv(2,prm2)-core->inv(2.,prm))/1e-8,              "linearCore inverse derivative(1)");
	prm2[1] = prm[1];
	delete core;

	x = new std::vector<double> (6,0);
	k = new std::vector<int> (6,0);
	n = new std::vector<int> (6,50);
	for (i=0; i<6; i++) (*x)[i] = 2*i+.1;
	(*k)[0] = 24; (*k)[1] = 32; (*k)[2] = 40; (*k)[3] = 48; (*k)[4] = 50; (*k)[5] = 48;
	data = new PsiData ( *x, *n, *k, 2 );
	core = new logCore ( data );
	th = exp(-2./3);

	failures += T->isequal ( core->g(th,prm), 0,                  "logCore at threshold");
	failures += T->isequal ( core->dgx(3.,prm),prm[0]/3.,        "linearCore derivative stimulus");
	failures += T->isequal ( core->dg(th,prm,0), log(th),         "logCore derivative(0) at threshold");
	failures += T->isequal ( core->dg(th,prm,1), 1,               "logCore derivative(1) at threshold");
	failures += T->isequal ( core->ddg(th,prm,0,0), 0,            "logCore derivative(0,0) at threshold");
	failures += T->isequal ( core->ddg(th,prm,1,0), 0,            "logCore derivative(1,0) at threshold");
	failures += T->isequal ( core->ddg(th,prm,1,1), 0,            "logCore derivative(1,1) at threshold");
	failures += T->isequal ( core->ddg(th,prm,0,1), core->ddg(th,prm,1,0), "logCore 2nd derivative symmetry at threshold");
	failures += T->isequal ( core->inv(0,prm), th,                "logCore inverse");
	failures += T->isequal ( core->inv(core->g(2.,prm),prm),2.,   "logCore inversion g(inv(2))");
	failures += T->isequal ( core->g(core->inv(2.,prm),prm),2.,   "logCore inversion inv(g(2))");
	prm2[0] += 1e-8;
	failures += T->isequal ( core->dinv(2.,prm,0), (core->inv(2,prm2)-core->inv(2,prm))/1e-8,            "logCore inversion dinv(2,0)");
	prm2[0] = prm[0];
	prm2[1] += 1e-8;
	failures += T->isequal ( core->dinv(2.,prm,1), (core->inv(2,prm2)-core->inv(2,prm))/1e-8, "logCore inversion dinv(2,1)");
	prm2[1] = prm[1];
	delete core;

	core = new NakaRushton ( data );
	prm2[0] = 4.; prm2[1] = 2.; prm2[2] = .02;
	prm = prm2;
	th = prm2[0];
	failures += T->isequal ( core->g(th,prm), 0.5,                  "NakaRushton at threshold");
	failures += T->isequal ( core->dgx(2.,prm), 96./400.,           "NakaRushton derivative stimulus");
	// Derivatives are tested separately
	failures += T->isequal ( core->inv(0.5,prm), th,                "NakaRushton inverse");
	failures += T->isequal ( core->inv(core->g(2.,prm),prm),2.,     "NakaRushton inversion inv(g(2))");
	failures += T->isequal ( core->g(core->inv(.5,prm),prm),.5,     "NakaRushton inversion g(inv(0.5))");
	prm2[0] += 1e-8;
	failures += T->isequal ( core->dinv(.2,prm,0), (core->inv(.2,prm2)-core->inv(.2,prm))/1e-8, "NakaRushton inversion dinv(2,0)");
	prm2[0] = prm[0];
	prm2[1] += 1e-8;
	failures += T->isequal ( core->dinv(.2,prm,1), (core->inv(.2,prm2)-core->inv(.2,prm))/1e-8, "NakaRushton inversion dinv(2,1)");
	prm2[1] = prm[1];
	delete core;


	delete data;
	delete x;
	delete k;
	delete n;

	return failures;
}

int LinalgTests ( TestSuite * T ) {
	// These tests compare the results with the respective numpy/scipy routines
	int failures (0);

	Matrix *M = new Matrix (3,3);
	std::vector<double> x(3),b(3);

	(*M)(0,0) = 0.75; (*M)(0,1) = 0.52; (*M)(0,2) = -.16;
	(*M)(1,0) = 0.52; (*M)(1,1) = 1.38; (*M)(1,2) = -.42;
	(*M)(2,0) = -.16; (*M)(2,1) = -.42; (*M)(2,2) = 0.53;

	Matrix *I = M->inverse();
	failures += T->isequal ( (*I)(0,0),  1.80488979, "Inverse (0,0)" );
	failures += T->isequal ( (*I)(1,0), -0.67772799, "Inverse (1,0)" );
	failures += T->isequal ( (*I)(2,0),  0.00780493, "Inverse (2,0)" );
	failures += T->isequal ( (*I)(0,1), -0.67772799, "Inverse (0,1)" );
	failures += T->isequal ( (*I)(1,1),  1.20943876, "Inverse (1,1)" );
	failures += T->isequal ( (*I)(2,1),  0.75382604, "Inverse (2,1)" );
	failures += T->isequal ( (*I)(0,2),  0.00780493, "Inverse (0,2)" );
	failures += T->isequal ( (*I)(1,2),  0.75382604, "Inverse (1,2)" );
	failures += T->isequal ( (*I)(2,2),  2.48652024, "Inverse (2,2)" );
	delete I;

	b[0] = 1; b[1] = 0.5; b[2] = 0;
	x = M->solve(b);
	failures += T->isequal ( x[0],  1.4660258,  "solving Ax=b, x[0]" );
	failures += T->isequal ( x[1], -0.0730086,  "solving Ax=b, x[1]" );
	failures += T->isequal ( x[2],  0.38471795, "solving Ax=b, x[2]" );

	I = M->cholesky_dec ();
	failures += T->isequal ( (*I)(0,0),  0.8660254,  "Cholesky (0,0)" );
	failures += T->isequal ( (*I)(1,0),  0.60044428, "Cholesky (1,0)" );
	failures += T->isequal ( (*I)(2,0), -0.18475209, "Cholesky (2,0)" );
	failures += T->isequal ( (*I)(0,1),  0.        , "Cholesky (0,1)" );
	failures += T->isequal ( (*I)(1,1),  1.00968642, "Cholesky (1,1)" );
	failures += T->isequal ( (*I)(2,1), -0.30610164, "Cholesky (2,1)" );
	failures += T->isequal ( (*I)(0,2),  0.        , "Cholesky (0,2)" );
	failures += T->isequal ( (*I)(1,2),  0.        , "Cholesky (1,2)" );
	failures += T->isequal ( (*I)(2,2),  0.63416753, "Cholesky (2,2)" );
	delete I;

	I = M->lu_dec ();
	failures += T->isequal ( (*I)(0,0),  0.75      , "LU (0,0)" );
	failures += T->isequal ( (*I)(1,0),  0.69333333, "LU (1,0)" );
	failures += T->isequal ( (*I)(2,0), -0.21333333, "LU (2,0)" );
	failures += T->isequal ( (*I)(0,1),  0.52      , "LU (0,1)" );
	failures += T->isequal ( (*I)(1,1),  1.01946667, "LU (1,1)" );
	failures += T->isequal ( (*I)(2,1), -0.30316505, "LU (2,1)" );
	failures += T->isequal ( (*I)(0,2), -0.16      , "LU (0,2)" );
	failures += T->isequal ( (*I)(1,2), -0.30906667, "LU (1,2)" );
	failures += T->isequal ( (*I)(2,2),  0.40216845, "LU (2,2)" );
	delete I;

	b = (*M)*x;
	failures += T->isequal ( b[0], 1., "Ax=b, b[0]" );
	failures += T->isequal ( b[1], .5, "Ax=b, b[1]" );
	failures += T->isequal ( b[2], 0., "Ax=b, b[2]" );

	T->isequal ( M->symmetric(), 1., "M should be symmetric" );

	int i,j;
	I = new Matrix (3,3);
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
		(*I)(i,j) = (*M)(i,j);
	I->scale(2);
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			failures += T->isequal ( (*I)(i,j), 2*(*M)(i,j), "matrix scaling" );
	delete I;

	// Test that should only give the right solution with pivoting
	(*M)(0,0) = 11; (*M)(0,1) = 44; (*M)(0,2) = 1;
	(*M)(1,0) = .1; (*M)(1,1) = .4; (*M)(1,2) = 3;
	(*M)(2,0) =  0; (*M)(2,1) =  1; (*M)(2,2) =-1;
	b[0] = b[1] = b[2] = 1;
	x = M->solve(b);
	// We need pivoting only to make sure that we are not blown off completely.
	failures += T->isequal ( x[0], -5.26445,  "pivot Ax=b, x[0]", .025 );
	failures += T->isequal ( x[1],  1.33131,  "pivot Ax=b, x[1]", .02 );
	failures += T->isequal ( x[2],  0.331307, "pivot Ax=b, x[2]", .02 );

	delete M;

	return failures;
}

int ReturnTest ( TestSuite * T ) {
	// In some cases, jackkifing doesn't terminate
	PsiCore * core = new mwCore ( NULL, 1, 0.1 );
	PsiSigmoid * sigmoid = new PsiLogistic ();
	PsiPsychometric *pmf = new PsiPsychometric ( 2, core, sigmoid );
	// std::vector<double> x (5);      x[0] = 1; x[1] = 2; x[2] = 3; x[3] = 4; x[4] = 5;
	// std::vector<int>    k (5,10);   k[1] = 9; k[2] = 8;
	// std::vector<int>    n (5,10);
	std::vector<double> x (3);      x[0] = 3.8091348774813367; x[1] = 4.3712635982077179; x[2] = 6.0913291693220737;
	std::vector<int>    k (3,20);   k[0] = 16;
	std::vector<int>    n (3,20);
	PsiData * data = new PsiData ( x, n, k, 2 );

	pid_t childpid;
	int hang;
	pid_t phang;

	if ((childpid = fork()) == 0) {
		// Run the optimizer in a child process with a limited amount of time
		PsiOptimizer *opt = new PsiOptimizer ( pmf, data );
		std::vector<double> solution ( opt->optimize ( pmf, data ) );
		delete opt;

		delete data;
		delete core;
		delete sigmoid;
		delete pmf;

		exit(0); // If we got this far, we kill the child process
 	} else {
		sleep(3); // We wait 3s for the child. Otherweise, we consider it as a failure.
		phang = waitpid ( childpid, &hang, WNOHANG );

		if (phang!=0) kill ( childpid, 9 ); // The optimizer process can be killed

		delete data;
		delete core;
		delete sigmoid;
		delete pmf;

		return T->isequal ( !phang, 0, "optimizer hung" );
	}
}

int GetstartTest ( TestSuite * T ) {
	int failures (0);

	/************************
	 * 2AFC
	 */
	std::vector <double> x ( 6 );
	std::vector <int>    n ( 6, 50 );
	std::vector <int>    k ( 6 );
	double xmin,xmax;
	double ymin,ymax;

	// Set up data
	x[0] =  0.; x[1] =  2.; x[2] =  4.; x[3] =  6.; x[4] =  8.; x[5] = 10.;
	k[0] = 3;  k[1] = 10;  k[2] = 34;  k[3] = 45;  k[4] = 50;  k[5] = 50;
	PsiData *data = new PsiData (x,n,k,1);

	// Set up psychometric function
	PsiPrior *prior = new UniformPrior(0.,0.1);
	abCore * core = new abCore();
	PsiLogistic * sigmoid = new PsiLogistic();
	PsiPsychometric *pmf = new PsiPsychometric ( 1, core, sigmoid );
	std::vector<double> prm(4);
	prm[0] = 4; prm[1] = 0.8; prm[2] = 0.02; prm[3] = 0.1;
	pmf->setPrior( 2, prior);
	pmf->setPrior( 3, prior);
	std::vector<double> start;
	start = getstart ( pmf, data, 7, 3, 3 );

	failures += T->isequal ( start[0], 3.33333,    "yes-no: Starting value for alpha", 1e-5 );
	failures += T->isequal ( start[1], 1.03704,   "yes-no: Starting value for beta", 1e-5 );
	failures += T->isequal ( start[2], 0.00185185, "yes-no: Starting value for lambda", 1e-5 );
	failures += T->isequal ( start[3], 0.0166667, "yes-no: Starting value for gamma", 1e-5 );

	a_range ( data, &xmin, &xmax );
	parameter_range ( data, pmf, 0, &ymin, &ymax );
	failures += T->isequal ( xmin,  0, "yes-no: minimum of alpha range", 1e-5 );
	failures += T->isequal ( xmax, 10, "yes-no: maximum of alpha range", 1e-5 );
	failures += T->isequal ( ymin,  0, "yes-no: minimum of alpha range", 1e-5 );
	failures += T->isequal ( ymax, 10, "yes-no: maximum of alpha range", 1e-5 );
	b_range ( data, &xmin, &xmax );
	parameter_range ( data, pmf, 1, &ymin, &ymax );
	failures += T->isequal ( xmin, 2.0, "yes-no: minimum of beta range", 1e-5 );
	failures += T->isequal ( xmax, 10.0, "yes-no: maximum of beta range", 1e-5 );
	failures += T->isequal ( ymin, 2.0, "yes-no: minimum of beta range", 1e-5 );
	failures += T->isequal ( ymax, 10.0, "yes-no: maximum of beta range", 1e-5 );
	lm_range ( data, &xmin, &xmax );
	parameter_range ( data, pmf, 2, &ymin, &ymax );
	failures += T->isequal ( xmin,  0, "yes-no: minimum of lambda range", 1e-5 );
	failures += T->isequal ( xmax, .1, "yes-no: maximum of lambda range", 1e-5 );
	failures += T->isequal ( ymin,  0, "yes-no: minimum of lambda range", 1e-5 );
	failures += T->isequal ( ymax, .1, "yes-no: maximum of lambda range", 1e-5 );
	gm_range ( data, &xmin, &xmax );
	parameter_range ( data, pmf, 3, &ymin, &ymax );
	failures += T->isequal ( xmin,  0, "yes-no: minimum of gamma range", 1e-5 );
	failures += T->isequal ( xmax, .1, "yes-no: maximum of gamma range", 1e-5 );
	failures += T->isequal ( ymin,  0, "yes-no: minimum of gamma range", 1e-5 );
	failures += T->isequal ( ymax, .1, "yes-no: maximum of gamma range", 1e-5 );

	delete data;
	delete pmf;

	// Set up data
	x[0] =  0.; x[1] =  2.; x[2] =  4.; x[3] =  6.; x[4] =  8.; x[5] = 10.;
	k[0] = 24;  k[1] = 32;  k[2] = 40;  k[3] = 48;  k[4] = 50;  k[5] = 48;
	data = new PsiData (x,n,k,2);

	// Set up psychometric function
	pmf = new PsiPsychometric ( 2, core, sigmoid );
	prm[0] = 4; prm[1] = 0.8; prm[2] = 0.02;
	pmf->setPrior( 2, prior );
	start = getstart( pmf, data, 7, 3, 3);

	failures += T->isequal ( start[0], 3.33333,  "2afc: Starting value for alpha", 1e-5 );
	failures += T->isequal ( start[1], 0.888889, "2afc: Starting value for beta", 1e-5 );
	failures += T->isequal ( start[2], 0.0203704, "2afc: Starting value for lambda", 1e-5 );

	a_range ( data, &xmin, &xmax );
	parameter_range ( data, pmf, 0, &ymin, &ymax );
	failures += T->isequal ( xmin,  0, "2afc: minimum of alpha range", 1e-5 );
	failures += T->isequal ( xmax, 10, "2afc: maximum of alpha range", 1e-5 );
	failures += T->isequal ( ymin,  0, "2afc: minimum of alpha range", 1e-5 );
	failures += T->isequal ( ymax, 10, "2afc: maximum of alpha range", 1e-5 );
	b_range ( data, &xmin, &xmax );
	parameter_range ( data, pmf, 1, &ymin, &ymax );
	failures += T->isequal ( xmin, 2.0, "2afc: minimum of beta range", 1e-5 );
	failures += T->isequal ( xmax, 10.0, "2afc: maximum of beta range", 1e-5 );
	failures += T->isequal ( ymin, 2.0, "2afc: minimum of beta range", 1e-5 );
	failures += T->isequal ( ymax, 10.0, "2afc: maximum of beta range", 1e-5 );
	lm_range ( data, &xmin, &xmax );
	parameter_range ( data, pmf, 2, &ymin, &ymax );
	failures += T->isequal ( xmin,  0, "2afc: minimum of lambda range", 1e-5 );
	failures += T->isequal ( xmax, .1, "2afc: maximum of lambda range", 1e-5 );
	failures += T->isequal ( ymin,  0, "2afc: minimum of lambda range", 1e-5 );
	failures += T->isequal ( ymax, .1, "2afc: maximum of lambda range", 1e-5 );

	delete core;
	delete sigmoid;
	delete prior;

	std::vector<double> pmin ( 3,  0 );
	std::vector<double> pmax ( 3, .5 );
	pmax[2] = .05;
	std::vector<double> u;
	char txt[200];
	unsigned int i;
	PsiGrid grid ( pmin, pmax, 5 );
	for ( i=0; i<3; i++ ) {
		sprintf ( txt, "grid parameter %d lower limit", i );
		failures += T->isequal ( grid.get_lower(i), pmin[i], txt, 1e-5 );
		sprintf ( txt, "grid parameter %d upper limit", i );
		failures += T->isequal ( grid.get_upper(i), pmax[i], txt, 1e-5 );
	}
	failures += T->isequal ( grid.empty(), false, "grid.empty() on nonempty grid" );
	failures += T->isequal ( PsiGrid().empty(), true, "grid.empty() on empty grid" );
	failures += T->isequal ( grid.dimension(), 3, "grid.dimension() on 3d grid" );
	failures += T->isequal ( grid.get_gridsize(), 5, "grid.get_gridsize() on small grid" );
	u = grid.front();
	failures += T->isequal ( u[0], 0.000, "grid.front()[0]" );
	failures += T->isequal ( u[1], 0.125, "grid.front()[1]" );
	failures += T->isequal ( u[2], 0.250, "grid.front()[2]" );
	failures += T->isequal ( u[3], 0.375, "grid.front()[3]" );
	failures += T->isequal ( u[4], 0.500, "grid.front()[4]" );

	PsiGrid newgrid;
	u = std::vector<double> ( 3, .125 );
	newgrid = grid.shift ( u );
	u = newgrid.front();
	failures += T->isequal ( u[0], -0.125,  "shifted grid.front()[0]" );
	failures += T->isequal ( u[1], 0.,      "shifted grid.front()[1]", 1e-5 );
	failures += T->isequal ( u[2], 0.125,   "shifted grid.front()[2]", 1e-5 );
	failures += T->isequal ( u[3], 0.25,    "shifted grid.front()[3]", 1e-5 );
	failures += T->isequal ( u[4], 0.375,   "shifted grid.front()[4]", 1e-5 );

	u = std::vector<double> ( 3, .125 );
	newgrid = grid.shrink ( u );
	u = newgrid.front ();
	failures += T->isequal ( u[0], 0.,      "shrunken grid.front()[0]" );
	failures += T->isequal ( u[1], 0.0625,  "shrunken grid.front()[1]", 1e-5 );
	failures += T->isequal ( u[2], 0.125,   "shrunken grid.front()[2]", 1e-5 );
	failures += T->isequal ( u[3], 0.1875,  "shrunken grid.front()[3]", 1e-5 );
	failures += T->isequal ( u[4], 0.25,    "shrunken grid.front()[4]", 1e-5 );

	newgrid = grid.subgrid ();
	failures += T->isequal ( newgrid.dimension(), 2, "subgrid dimension" );

	u = linspace ( 0,1,3 );
	failures += T->isequal ( u[0], 0,   "linspace 0" );
	failures += T->isequal ( u[1], 0.5, "linspace 1" );
	failures += T->isequal ( u[2], 1,   "linspace 2" );

	std::list< std::vector<double> > gridpoints;
	std::list< std::vector<double> >::iterator i_gp;
	grid = PsiGrid ( pmin, pmax, 2 );
	makegridpoints ( grid, u, 0, &gridpoints );
	failures += T->isequal ( gridpoints.size(), 8, "Number of generated gridpoints from 2x2x2 grid" );
	for ( i=0, i_gp=gridpoints.begin(); i_gp!=gridpoints.end(); i_gp++, i++ ) {
		sprintf ( txt, "gridpoint %d first param", i );
		failures += T->isequal ( (*i_gp)[0], (i<4 ? 0 : 0.5), txt );
		sprintf ( txt, "gridpoint %d second param", i );
		failures += T->isequal ( (*i_gp)[1], ((i/2)%2 == 0 ? 0 : 0.5), txt );
		sprintf ( txt, "gridpoint %d third param", i );
		failures += T->isequal ( (*i_gp)[2], (i%2 == 0 ? 0 : 0.05), txt );
	}

	std::list< std::vector<double> > bestprm;
	std::list< double > L;
	evalgridpoints ( gridpoints, &bestprm, &L, data, pmf, 2 );


	failures += T->isequal ( L.front(), 35.4381, "Best fit on grid", 1e-4 );
	failures += T->isequal ( L.back(),  43.346 , "Second best fit on grid", 1e-4 );
	failures += T->isequal ( bestprm.front()[0], .5,  "Best fitting alpha on grid" );
	failures += T->isequal ( bestprm.front()[1], .5,  "Best fitting beta on grid" );
	failures += T->isequal ( bestprm.front()[2], .05, "Best fitting lambda on grid" );
	failures += T->isequal ( bestprm.back()[0], .0,   "Second best fitting alpha on grid" );
	failures += T->isequal ( bestprm.back()[1], .5,   "Second best fitting beta on grid" );
	failures += T->isequal ( bestprm.back()[2], .05,  "Second best fitting lambda on grid" );

	gridpoints = std::list< std::vector<double> > (0);
	std::list< PsiGrid > newgrids;
	bestprm.pop_back();  // Delete the last element to keep number of points small
	updategridpoints ( grid, bestprm, &gridpoints, &newgrids );
	for ( i=0, i_gp=gridpoints.begin(); i_gp!=gridpoints.end(); i_gp++, i++ ) {
		sprintf ( txt, "gridpoint %d first param", i );
		failures += T->isequal ( (*i_gp)[0], (i<4 ? 0.25 : 0.75), txt );
		sprintf ( txt, "gridpoint %d second param", i );
		failures += T->isequal ( (*i_gp)[1], ((i/2)%2 == 0 ? .25 : 0.75), txt );
		sprintf ( txt, "gridpoint %d third param", i );
		failures += T->isequal ( (*i_gp)[2], (i%2 == 0 ? .025 : 0.075), txt );
	}

	delete data;
	delete pmf;

	return failures;
}

int IntegrateTest ( TestSuite * T ) {
	int failures (0);
	unsigned int i;

	std::vector<double> intensity ( 6 );
	std::vector<int>    ncorrect  ( 6 );
	std::vector<int>    ntrials   ( 6, 50 );

	intensity[0] =  0.;    ncorrect[0]  = 24;
	intensity[1] =  2.;    ncorrect[1]  = 32;
	intensity[2] =  4.;    ncorrect[2]  = 40;
	intensity[3] =  6.;    ncorrect[3]  = 48;
	intensity[4] =  8.;    ncorrect[4]  = 50;
	intensity[5] = 10.;    ncorrect[5]  = 48;
	PsiData * data = new PsiData ( intensity, ntrials, ncorrect, 2 );

	mwCore * core = new mwCore();
	PsiLogistic * sigmoid = new PsiLogistic();
	PsiPsychometric * pmf = new PsiPsychometric ( 2, core, sigmoid );

	//std::vector<double> grid;
	//std::vector<double> tar ( 7 );
	//std::vector<double> start ( 3 );
	//grid = raw_grid ( data, pmf, 0, 7 );
	//failures += T->isequal ( grid[0], 0, "Lower bound of m" );
	//failures += T->isequal ( grid[6], 10, "Upper bound of m" );
	//failures += T->isequal ( grid.size(), 7, "Length of m grid" );
	//grid = raw_grid ( data, pmf, 1, 7 );
	//failures += T->isequal ( grid[0], 2, "Lower bound of w" );
	//failures += T->isequal ( grid[6], 10, "Upper bound of w" );
	//failures += T->isequal ( grid.size(), 7, "Length of w grid" );
	//grid = raw_grid ( data, pmf, 2, 7 );
	//failures += T->isequal ( grid[0], 0, "Lower bound of lm" );
	//failures += T->isequal ( grid[6], 0.1, "Upper bound of lm" );
	//failures += T->isequal ( grid.size(), 7, "Length of lm grid" );

	//PsiPrior * dist = new GaussPrior ( 0, 1 );
	//tar[0] = -1.28155157;
	//tar[1] = -0.72791329;
	//tar[2] = -0.34069483;
	//tar[3] = 0.;
	//tar[4] =  0.34069483;
	//tar[5] = 0.72791329;
	//tar[6] =  1.28155157;
	//grid = cdf_grid ( dist, 0.1, 0.9, 7 );
	//for ( i=0; i<7; i++ ) {
	//	failures += T->isequal ( grid[i], tar[i], "cdf_grid Gauss", 1e-4 );
	//	tar[i] = dist->pdf ( grid[i] );
	//}
	//start[0] = -2; start[1] = 1.5; start[2] = -9;
	//start = fit_posterior ( grid, tar, start, 0 );
	//failures += T->isequal ( start[0], 0, "fit_posterior Gauss mean", 1e-4 );
	//failures += T->isequal ( start[1], 1, "fit_posterior Gauss std", 1e-4 );
	//delete dist;

	//dist = new GammaPrior ( 4, 1 );
	//tar[0] = 1.74476956;
	//tar[1] = 2.4572856;
	//tar[2] = 3.06204868;
	//tar[3] = 3.67206075;
	//tar[4] = 4.35885352;
	//tar[5] = 5.23689165;
	//tar[6] = 6.68078307;
	//grid = cdf_grid ( dist, 0.1, 0.9, 7 );
	//for ( i=0; i<7; i++ ) {
	//	failures += T->isequal ( grid[i], tar[i], "cdf_grid Gamma", 1e-4 );
	//	tar[i] = dist->pdf ( grid[i] );
	//}
	//start[0] = 2; start[1] = 1.5; start[2] = -9;
	//start = fit_posterior ( grid, tar, start, 1 );
	//failures += T->isequal ( start[0], 4, "fit_posterior Gamma shape", 1e-3 );
	//failures += T->isequal ( start[1], 1, "fit_posterior Gamma scale", 1e-3 );
	//delete dist;

	//dist = new BetaPrior ( 2, 20 );
	//tar[0] = 0.025617;
	//tar[1] = 0.043710;
	//tar[2] = 0.060631;
	//tar[3] = 0.078644;
	//tar[4] = 0.099659;
	//tar[5] = 0.127175;
	//tar[6] = 0.172935;
	//grid = cdf_grid ( dist, 0.1, 0.9, 7 );
	//for ( i=0; i<7; i++ ) {
	//	failures += T->isequal ( grid[i], tar[i], "cdf_grid Beta", 1e-4 );
	//	tar[i] = dist->pdf ( grid[i] );
	//}
	//start = fit_posterior ( grid, tar, start, 2 );
	//failures += T->isequal ( start[0], 2, "fit_posterior Beta alpha", 1e-3 );
	//failures += T->isequal ( start[1], 20, "fit_posterior Beta beta", 1e-2 );
	//delete dist;

	//PsiIndependentPosterior posterior = independent_marginals ( pmf, data, 3, 7 );

	//failures += T->isequal ( posterior.get_posterior ( 0 )->getprm ( 0 ), 3.280,  "Posterior for m -- mu", 1e-3 );
	//failures += T->isequal ( posterior.get_posterior ( 0 )->getprm ( 1 ), 0.436,  "Posterior for m -- sg", 1e-3 );
	//failures += T->isequal ( posterior.get_posterior ( 1 )->getprm ( 0 ), 12.732, "Posterior for w -- k", 1e-3 );
	//failures += T->isequal ( posterior.get_posterior ( 1 )->getprm ( 1 ), 0.362,  "Posterior for w -- th", 1e-3 );
	//failures += T->isequal ( posterior.get_posterior ( 2 )->getprm ( 0 ), 3.170,  "Posterior for lm -- al", 1e-3 );
	//failures += T->isequal ( posterior.get_posterior ( 2 )->getprm ( 1 ), 111.26, "Posterior for lm -- bt", 1e-3 );

	//MCMCList samples = sample_posterior ( pmf, data, posterior, 600 );
	//failures += T->isequal ( samples.getMean ( 0 ), posterior.get_posterior ( 0 )->getprm ( 0 ), "Sampled and fitted posterior mean for m", .2 );
	//failures += T->isequal ( samples.getMean ( 1 ),
	//		posterior.get_posterior ( 1 )->getprm ( 0 ) * posterior.get_posterior ( 1 )->getprm(1),
	//		"Sampled and fitted posterior mean for w", 1e-2 );
	//failures += T->isequal ( samples.getMean ( 2 ),
	//		posterior.get_posterior ( 2 )->getprm ( 0 ) / ( posterior.get_posterior ( 2 )->getprm(0)+posterior.get_posterior ( 2 )->getprm(1) ),
	//		"Sampled and fitted posterior mean for lm", 1e-2 );
	//failures += T->isless ( samples.get_accept_rate (), 0.1, "number of duplicates in SIR" );

	//sample_diagnostics ( pmf, data, &samples );
	//// Simply check whether these values have been set (i.e. whether they are not 0)
	//failures += T->isequal ( samples.getRkd ( 0 )!=0, true, "Rkd is set" );
	//failures += T->isequal ( samples.getRpd ( 0 )!=0, true, "Rpd is set" );
	//failures += T->isequal ( samples.getlogratio ( 0, 0 )!=0, true, "logratio is set" );
	//failures += T->isequal ( samples.getppRkd ( 0 )!=0, true, "ppRkd is set" );
	//failures += T->isequal ( samples.getppRpd ( 0 )!=0, true, "ppRpd is set" );
	//failures += T->isequal ( samples.getppData ( 0 ).size(), 6, "ppData is set" );
	//failures += T->isequal ( samples.getdeviance ( 0 )!= 0, true, "Deviance is set" );
	//failures += T->isequal ( samples.getppDeviance ( 0 )!= 0, true, "ppDeviance is set" );

	/*  This is if you want to get the fits to plot them with gnuplot
	double x;
	unsigned int myprm ( 0 );
	std::vector<double> wgrid;
	std::vector<double> wmargin;
	for ( myprm=0; myprm<3; myprm++ ) {
		wmargin = posterior.get_margin ( myprm );
		wgrid = posterior.get_grid ( myprm );
		x = posterior.get_posterior ( myprm )->pdf ( wgrid[4] ) / wmargin[4];
		for ( i=0; i<7; i++ ) {
			wmargin[i] *= x;
			std::cout << wgrid[i] << " " << wmargin[i] << "\n";
		}
		std::cout << "\n\n";

		if ( myprm < 2 )
			for ( x=0; x<10; x+=.1 ) {
				std::cout << x << " " << posterior.get_posterior ( myprm )->pdf ( x ) << "\n";
			}
		else
			for ( x=0; x<.1; x+=.001 ) {
				std::cout << x << " " << posterior.get_posterior ( myprm )->pdf ( x ) << "\n";
			}
		std::cout << "\n\n";
	}
	*/

	return failures;
}

int main ( int argc, char ** argv ) {
	TestSuite Tests ( "tests_all.log" );
	Tests.addTest(&PsychometricValues,    "Values of the psychometric function");
	// Tests.addTest(&BetaModelTest,         "Beta psychometric function model");
	Tests.addTest(&DerivativeCheck,       "Derivaties of elements" );
	Tests.addTest(&OptimizerSolution,     "Solutions of optimizer");
	Tests.addTest(&BootstrapTest,         "Bootstrap properties");
	Tests.addTest(&SigmoidTests,          "Properties of sigmoids");
	Tests.addTest(&CoreTests,             "Tests of core objects");
	Tests.addTest(&MCMCTest,              "MCMC");
	Tests.addTest(&PriorTest,             "Priors");
	Tests.addTest(&LinalgTests,           "Linear algebra routines");
	Tests.addTest(&ReturnTest,            "Testing return bug in jackknifedata");
	Tests.addTest(&InitialParametersTest, "Initial parameter heuristics" );
	Tests.addTest(&GetstartTest,          "Finding good starting values" );
	Tests.addTest ( &IntegrateTest,        "Approximate numerical integration" );

	Tests.runTests();
}
