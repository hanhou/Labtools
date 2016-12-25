/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#include "optimizer.h"
#include "getstart.h"
#include <cmath>
#include <limits>

// #define DEBUG_OPTIMIZER

#ifdef DEBUG_OPTIMIZER
#include <iostream>
#include <fstream>
std::ofstream logfile ( "optimizer.log" );
#endif

const double maxstep (1e-7);
const double maxfstep(1e-7);
const int    maxiter (80);

PsiOptimizer::PsiOptimizer ( const PsiPsychometric * model, const PsiData * data)
	: nparameters ( model->getNparams() ),
	simplex ( nparameters+1, std::vector<double> (nparameters) ),fx ( nparameters+1 ),
	x  ( nparameters ),
	xx ( nparameters ),
	start ( nparameters ),
	modified ( nparameters+1, true )
{}

PsiOptimizer::~PsiOptimizer ( void ) {}

double testfunction(const std::vector<double>& x) {
	double out(0);
	unsigned int k;
	for (k=0; k<x.size(); k++)
		out += x[k]*x[k];
	return out;
}

/* While fixing broken unit tests we discovered that the Psignifit3
 * optimizer sometimes has problems reaching the same optimum as the
 * Psignifit2 optimizer. I particular this happens when lambda is near zero,
 * since the simplex is very close to a region that gives an infinite value
 * for the error function. This can be imagined as though the optimizer has
 * trouble climbing down a straight wall imposed by the hard constraint that
 * lambda must be in the open interval (0, 1). The problem is slight difference
 * in deviance between the values.
 *
 * To circumvent this we use the following transformation of the lambda variable
 * during optimization:
 *
 * lambda_hat = log(lambda/1-lambda)              <-------- logit
 *
 * lambda = (1/(1+exp(-lambda_hat)                <-------- logistic
 *
 * This effectively maps the value of lambda to lambda_hat and back. Lambda is
 * in the open interval (0, 1), whereas lambda_hat is in the space of real
 * numbers. This should make it much easier for the simplex, since no
 * constraints are imposed on the value lambda_hat. Some initial testing shows
 * that the Psignifit3 optimizer now approaches the Psignifit2  solution
 * closer. The transformation lamda -> lambda_hat takes place at the beginning
 * of the optimization run. The transformation lambda_hat -> lambda happens
 * before each evaluation of the error function (negloglikelihood) and when
 * returning the final value.
 *
 * While we were here, we did the same for gamma, since it is a rate it should
 * also never be less than 0 or greater than 1.
 *
 */

double lgst ( double x ) {
	return 1./(1+exp(-x));
}
double lgit ( double p ) {
	return log ( p/(1-p) );
}

void copy_lgst(const std::vector<double>& in, std::vector<double>& out, int nparameters){
	int l;
	for ( l=0; l<nparameters; l++ ) {
		out[l] = in[l];
		if ( l==2 || l==3 ) {
			out[l] = lgst ( out[l] );
		}
	}
}


std::vector<double> PsiOptimizer::optimize ( const PsiPsychometric * model, const PsiData * data, const std::vector<double>* startingvalue )
{
	int k, l;
	std::vector<double> incr ( model->getNparams() );
	if (startingvalue==NULL) {
		// start = model->getStart(data);
		start = getstart ( model, data, 8, 3, 3, &incr );
	} else {
		start = std::vector<double>(model->getNparams());
		incr  = std::vector<double>(model->getNparams());
		for ( k=0; k<int(model->getNparams()); k++ ) {
			start[k] = startingvalue->at(k);
			if ( (k+model->getNparams())<startingvalue->size() ) {
				incr[k] = startingvalue->at(k+model->getNparams());
			} else {
				incr[k] = 0.1 * start[k];
			}
		}
	}

	for ( k=0; k<nparameters+1; k++ ) {
		for ( l=0; l<nparameters; l++)
			simplex[k][l] = start[l];
		modified[k] = true;
	}

#ifdef DEBUG_OPTIMIZER
	logfile << "Starting values for optimization:\n";
	for (k=0; k<nparameters; k++)
		logfile << start[k] << "\n";
	logfile.flush();
#endif


	double ffx;         // function value at potential new point
	int maxind(0);      // Index of simplex node with maximal function value
	int minind(0);      // Index of simplex node with minimal function value
	double stepsize(1); // Measure for the size of the step
	double fstepsize(1);// Measure of size of a step in function values
	int iter(0);        // Number of simplex iterations
	int run;            // the model should be rerun after convergence
	double d;
	std::vector<double> output ( start );
	std::vector<double> prm ( start );


	for (run=0; run<2; run++) {
		for (k=1; k<nparameters+1; k++) {
			d = incr[k-1];
			simplex[k][k-1] += d;
			if ( model->evalPrior ( k-1, simplex[k][k-1] ) > 1000 ) {
				simplex[k][k-1] -= 2*d;
			}
#ifdef DEBUG_OPTIMIZER
			std::cout << "Start (regular," << k << "): " << simplex[k][0] << " " << simplex[k][1] << " " << simplex[k][2];
			if (nparameters>3)
				std::cout << " " << simplex[k][3];
			std::cout << "\n";
#endif
		}



		// transform starting values to logit
		for ( k=0; k<nparameters+1; k++ ) {
			simplex[k][2] = lgit ( simplex[k][2] );
			if ( nparameters > 3 ) {
				simplex[k][3] = lgit ( simplex[k][3] );
			}
		}

#ifdef DEBUG_OPTIMIZER
		std::cout << "Start (logit): " << simplex[0][0] << " " << simplex[0][1] << " " << simplex[0][2];
		if (nparameters>3)
			std::cout << " " << simplex[0][3];
		std::cout << "\n";
#endif

		// for (k=1; k<nparameters+1; k++) simplex[k][k-1] += .05;
		iter = 0;
		while (1) {
			// Evaluate model at every simplex node and determine maximum and minimum
			maxind = minind = 0;
			for (k=0; k<nparameters+1; k++) {
				if (modified[k]) {
					copy_lgst(simplex[k], prm, nparameters);
					fx[k] = model->neglpost(prm, data );
					modified[k] = false;
				}
				// fx[k] = testfunction(simplex[k]);
#ifdef DEBUG_OPTIMIZER
				for (l=0; l<nparameters; l++)
					logfile << simplex[k][l] << " ";
				logfile << fx[k] << "\n";
				logfile.flush();
#endif
				if (fx[k]<fx[minind]) minind = k;
				if (fx[k]>fx[maxind]) maxind = k;
			}

			// Avoid inf
			for ( k=0; k<nparameters+1; k++ ) {
				if ( fx[k] == std::numeric_limits<double>::infinity() ) {
					for ( l=0; l<nparameters; l++ ) {
						simplex[k][l] = start[l];
					}
					copy_lgst(simplex[k], prm, nparameters);
					fx[k] = model->neglpost(prm, data );
				}
			}

			// Check Stoping criteria based on simplex and function values
			stepsize = 0;
			for (k=0; k<nparameters; k++)
				stepsize += (simplex[maxind][k]-simplex[minind][k])*(simplex[maxind][k]-simplex[minind][k]);
			// Simplex size
			if (stepsize<maxstep) {
#ifdef DEBUG_OPTIMIZER
				logfile << "Terminating optimization due to small simplex size (" << stepsize << ") after " << iter << " iterations\n";
				logfile.flush();
#endif
				break;
			}
			// function value differences
			if ((fstepsize=(fx[maxind]-fx[minind])) < maxfstep ) {
#ifdef DEBUG_OPTIMIZER
				logfile << "Terminating optimization due to small function value variation (" << fstepsize << ") after " << iter << " iterations\n";
				logfile.flush();
#endif
				break;
			}

#ifdef DEBUG_OPTIMIZER
			logfile << iter << " " << fx[minind] << " " << stepsize << "\n";
			logfile.flush();
#endif

			// Calculate the average of the non maximal nodes
			for (k=0; k<nparameters; k++) x[k] = 0;
			for (k=0; k<nparameters+1; k++) {
				if (k!=maxind)
					for (l=0; l<nparameters; l++)
						x[l] += simplex[k][l];
			}
			for (k=0; k<nparameters; k++) x[k] /= nparameters;

			// Determine the reflection of the worst point
			for (k=0; k<nparameters; k++) xx[k] = x[k] - (simplex[maxind][k]-x[k]);

			// Now check what to do
			copy_lgst(xx, prm, nparameters);
			ffx = model->neglpost(prm,data);
			// ffx = testfunction(xx);
			if (ffx<fx[minind]) {
				// The reflected point is better than the previous worst point ~> Expand
				for (k=0; k<nparameters; k++) simplex[maxind][k] = x[k] - 2*(simplex[maxind][k] - x[k]);
				modified[maxind] = true;
			} else if (ffx>fx[maxind]) {
				// The reflected point is even worse than it was before ~> Shrink
				for (k=0; k<nparameters+1; k++) {
					for (l=0; l<nparameters; l++)
						simplex[k][l] = simplex[minind][l] + 0.5 * (simplex[k][l] - simplex[minind][l]);
					modified[k] = true;
				}
			} else {
				// The reflected point is somewhere in between
				for (k=0; k<nparameters; k++) simplex[maxind][k] = xx[k];
				fx[maxind] = ffx;
			}

			// Also cancel if the number of iterations gets to large
			if (iter++ > maxiter) {
#ifdef DEBUG_OPTIMIZER
				logfile << "Terminating optimization due to large number of iterations (" << iter << "). Final stepsize: " << stepsize << "\n";
				logfile.flush();
#endif
				break;
			}
		}

		// Evaluate model at every simplex node and determine minimum
		minind = 0;
		for (k=0; k<nparameters+1; k++) {
			if (modified[k]) {
				for ( l=0; l<nparameters; l++ ) {
					prm[l] = simplex[k][l];
					if ( l==2 || l==3 ) {
						prm[l] = lgst ( prm[l] );
					}
				}
				fx[k] = model->neglpost( prm, data );
				modified[k] = false;
			}
			// fx[k] = testfunction(simplex[k]);
			if (fx[k]<fx[minind]) minind = k;
		}

#ifdef DEBUG_OPTIMIZER
		std::cout << "Stop (logit): " << simplex[minind][0] << " " << simplex[minind][1] << " " << simplex[minind][2];
		if (nparameters>3)
			std::cout << " " << simplex[minind][3];
		std::cout << "\n";
#endif

		for ( k=0; k<nparameters+1; k++ ) {
			simplex[k][2] = lgst ( simplex[k][2] );
			if ( nparameters>3 ) {
				simplex[k][3] = lgst ( simplex[k][3] );
			}
		}

#ifdef DEBUG_OPTIMIZER
		std::cout << "Stop (regular): " << simplex[minind][0] << " " << simplex[minind][1] << " " << simplex[minind][2];
		if (nparameters>3)
			std::cout << " " << simplex[minind][3];
		std::cout << "\n";
#endif

		for (l=0; l<nparameters; l++) {
			output[l] = simplex[minind][l];
			// output[k] = start[k];
			// simplex[minind][l] = simplex[minind][0];
			simplex[0][l] = output[l];
			modified[l] = true;
		}
	}

	/*
	// Perform some Gradient descent steps
	for (k=0; k<40; k++) {
		x = model->dnegllikeli ( start, data );
		dl = 0;
		for (l=0; l<nparameters; l++) {
			start[l] -= .1*x[l];
			if (fabs(x[l])>dl)
				dl = fabs(x[l]);
		}
		if (dl<1e-6)
			break;
	}
	*/

#ifdef DEBUG_OPTIMIZER
	logfile << "Returning\n"; logfile.flush();

	for (l=0; l<data->getNblocks(); l++)
		logfile << data->getIntensity(l) << " " << data->getNcorrect(l) << " " << data->getPcorrect(l) << " " << data->getNtrials(l) << "\n";

	for (l=0; l<nparameters; l++) logfile << output[l] << " ";
	logfile.flush();
	logfile << "\n"; logfile.flush();
	logfile << "Done\n"; logfile.flush();

#endif

	return output;
}
