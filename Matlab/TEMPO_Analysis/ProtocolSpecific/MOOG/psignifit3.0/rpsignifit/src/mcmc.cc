/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#include "mcmc.h"

// #define DEBUG_MCMC

#include <iostream>
#include <iomanip>

/**********************************************************************
 *
 * MetropolisHastings sampling
 *
 */

MetropolisHastings::MetropolisHastings ( const PsiPsychometric * pmf, const PsiData * dat, PsiRandom * proposal )
	: PsiSampler ( pmf, dat ),
	propose(proposal->clone()),
	currenttheta(pmf->getNparams(),0),
	newtheta(pmf->getNparams(),0),
	stepwidths(pmf->getNparams(),.1),
	accept(0),
	qold(-1e5)
{
#ifdef DEBUG_MCMC
	std::cerr << "Hi my name is MetropolisHastings\n";
#endif
	setTheta ( currenttheta );
    currentdeviance = (pmf->deviance(currenttheta,dat));
}

std::vector<double> MetropolisHastings::draw ( void ) {
	double qnew, acc(propose->rngcall());
	const PsiPsychometric * model (getModel());
	const PsiData * data (getData());

	// propose a new point
	proposePoint(currenttheta, stepwidths, propose, newtheta);

	// negative log posterior of the point
	qnew = acceptance_probability ( currenttheta, newtheta );
	// std::cerr << qnew-qold << " " << exp(qnew-qold) << "\n";

	if (acc<exp(qnew-qold)) {
		// accept the new point
		qold = qnew;
		currenttheta = newtheta;
		currentdeviance = model->deviance ( currenttheta, data );
		accept ++;
#ifdef DEBUG_MCMC
		std::cerr << " ACCEPTED ";
#endif
	}
#ifdef DEBUG_MCMC
	else
		std::cerr << " REJECTED ";


	std::cout << "\n";
#endif

	return currenttheta;
}

double MetropolisHastings::acceptance_probability (
		const std::vector<double>& current_theta,
		const std::vector<double>& new_theta ) {
	double qnew;

	qnew = -getModel()->neglpost ( new_theta, getData() );

#ifdef DEBUG_MCMC

	std::cerr
		<< "Q_old: " << std::setiosflags ( std::ios::fixed ) << qold
		<< " Q_new: " << std::setiosflags ( std::ios::fixed ) << qnew
		<< " P(accept):" << std::setiosflags ( std::ios::fixed ) << (lratio>1 ? 1 : lratio)
        << "\n";
    int i, Nparams(getModel()->getNparams());
    std::cerr << "Current Theta:\t\t";
    for(i=0; i<Nparams; i++){
        std::cerr << currenttheta[i] << "\t";
    }
    std::cerr << "\n";

    std::cerr << "New Theta:\t\t";
    for(i=0; i<Nparams; i++){
        std::cerr << newtheta[i] << "\t";
    }
    std::cerr << "\n";

#endif

	return qnew;
}

void MetropolisHastings::proposePoint( std::vector<double> &current_theta,
										std::vector<double> &step_widths,
										PsiRandom * proposal,
										std::vector<double> &new_theta){
	const PsiPsychometric * model( getModel() );
	int prm, Nprm(model->getNparams());
	for (prm=0; prm<Nprm; prm++) {
		new_theta[prm] = current_theta[prm] + step_widths[prm] * proposal->draw ( );
	}
}

void MetropolisHastings::setTheta ( const std::vector<double>& prm ) {
	if (prm.size()==currenttheta.size())
		currenttheta = prm;
	else
		throw BadArgumentError();
	qold = getModel()->neglpost( currenttheta, getData() );
}

void MetropolisHastings::setStepSize ( double size, unsigned int param ) {
	if ( param<getModel()->getNparams() )
		stepwidths[param] = size;
	else
		throw BadIndexError();
}

void MetropolisHastings::setStepSize ( const std::vector<double>& sizes ) {
	unsigned int i;
	for (i=0; i<stepwidths.size(); i++)
		stepwidths[i] = sizes[i];
}

MCMCList MetropolisHastings::sample ( unsigned int N ) {
	const PsiData * data ( getData() );
	const PsiPsychometric * model ( getModel() );
	accept = 0;
	MCMCList out ( N, model->getNparams(), data->getNblocks() );
	PsiData *localdata = new PsiData ( data->getIntensities(), data->getNtrials(), data->getNcorrect(), data->getNalternatives() );
	std::vector< PsiData* > reduceddata (data->getNblocks() );
	std::vector<int> posterior_predictive ( data->getNblocks() );
	std::vector<double> probs ( data->getNblocks() );
	std::vector<double> est ( model->getNparams() );
	unsigned int i,j,k,l;

	std::vector<double> reducedx ( data->getNblocks()-1 );
	std::vector<int> reducedk ( data->getNblocks()-1 );
	std::vector<int> reducedn ( data->getNblocks()-1 );

	for ( k=0; k<data->getNblocks(); k++ ) {
		j = 0;
		for ( l=0; l<data->getNblocks(); l++ ) {
			if ( l!=k ) {
				reducedx[j] = data->getIntensity(l);
				reducedk[j] = data->getNcorrect(l);
				reducedn[j] = data->getNtrials(l);
				j++;
			}
		}
		reduceddata[k] = new PsiData ( reducedx, reducedn, reducedk, data->getNalternatives() );
	}

	qold = acceptance_probability ( currenttheta, currenttheta );

	for (i=0; i<N; i++) {
		// Draw the next sample
		est = draw();
		out.setEst ( i, est, 0. );
		out.setdeviance ( i, getDeviance() );

		// determine posterior predictives
		for ( k=0; k<data->getNblocks(); k++ )
			probs[k] = model->evaluate ( data->getIntensity(k), est );
		newsample ( localdata, probs, &posterior_predictive);
		localdata->setNcorrect ( posterior_predictive );
		out.setppData ( i, posterior_predictive, model->deviance ( est, localdata ) );

		probs = model->getDevianceResiduals ( est, data );
		out.setRpd ( i, model->getRpd ( probs, est, data ) );
		out.setRkd ( i, model->getRkd ( probs, data ) );

		probs = model->getDevianceResiduals ( est, localdata );
		out.setppRpd ( i, model->getRpd ( probs, est, localdata ) );
		out.setppRkd ( i, model->getRkd ( probs, localdata ) );

		// Store log posterior ratios for reduced data sets
		for ( k=0; k<data->getNblocks(); k++) {
			/*
			j=0;
			for ( l=0; l<data->getNblocks(); l++ ) {
				if ( l!=k ) {
					reducedx[j] = data->getIntensity(l);
					reducedk[j] = data->getNcorrect(l);
					reducedn[j] = data->getNtrials(l);
					j++;
				}
			}
			reduceddata = new PsiData ( reducedx, reducedn, reducedk, data->getNalternatives() );
			*/
			out.setlogratio ( i, k, model->neglpost(est,data)-model->neglpost(est,reduceddata[k]) );
		}
#ifdef DEBUG_MCMC
		std::cerr << " accept: " << std::setiosflags ( std::ios::fixed ) << double(accept)/(i+1) << "\n";
#endif
	}

#ifdef DEBUG_MCMC
	std::cerr << "Acceptance rate: " << double(accept)/N << "\n";
#endif
	out.set_accept_rate(double(accept)/N);

	delete localdata;
	for ( k=0; k<reduceddata.size(); k++ ) {
		delete reduceddata[k];
	}

	return out;
}


/**********************************************************************
 *
 * Generic Metropolis MCMC
 *
 */

void GenericMetropolis::proposePoint(std::vector<double> &current_theta,
									  std::vector<double> &step_widths,
									  PsiRandom * proposal,
									  std::vector<double> &new_theta) {
	const PsiPsychometric * model ( getModel() );

	/* update one direction of theta */
	new_theta = current_theta;
	new_theta[currentindex] += step_widths[currentindex] * proposal->draw();

	/* iterate parameter index each time a new point is proposed */
	currentindex = (currentindex + 1) % model->getNparams();
}


void GenericMetropolis::findOptimalStepwidth( PsiMClist const &pilot ){
    if ( pilot.getNsamples() < pilot.getNparams() +1 ){
        throw BadArgumentError("The number of samples in the pilot must be at least equal to the number of free parameters.");
    }
	int i,j,prm, Nparams(pilot.getNparams()), Nsamples(pilot.getNsamples());
	double std_residuals; // standard deviation of the residuals
	int *paramindex = new int[Nparams-1];
	Matrix X = Matrix(Nsamples, Nparams+1); // extended data matrix

	for (prm=0; prm<Nparams; prm++){
		for (i=0; i<prm; i++) paramindex[i] = i;
		for (i=prm+1; i<Nparams; i++) paramindex[i-1] = i;

		/* copy sampled data in a matrix: */
		for (i=0; i<Nsamples; i++){ // iterate samples
			X(i,0) = 1.0; // fill first column with 1.
			for (j=0; j<Nparams-1; j++){
				X(i,j+1) = pilot.getEst(i,paramindex[j]);
			}
			X(i,Nparams) = pilot.getEst(i,prm); // last column is the target
		}

		/* do QR-decomposition: */
		Matrix *R = X.qr_dec();

		/* get the residuals from R and calculate the STD: */
		std_residuals = sqrt( (*R)(Nparams,Nparams) * (*R)(Nparams,Nparams) / double(Nsamples) );

		/* multiply std deviation with 2.38/sqrt(Nparams) as suggested by Gelman et al. (1995) */
		setStepSize( std_residuals * 2.38 / sqrt(double(Nparams)), prm );

		delete R;
	}
	delete [] paramindex;
}

/**********************************************************************
 *
 * DefaultMCMC
 *
 */
DefaultMCMC::DefaultMCMC ( const PsiPsychometric* Model, const PsiData* Data, PsiRandom* prop ) :
	MetropolisHastings ( Model, Data, new GaussRandom ),
	proposaldistributions ( Model->getNparams () )
{
#ifdef DEBUG_MCMC
	std::cerr << "Hi my name is DefaultMCMC\n";
#endif
}

DefaultMCMC::~DefaultMCMC ( void ) {
	unsigned int i;
	for (i=0; i<proposaldistributions.size(); i++) {
		delete proposaldistributions[i];
	}
}

double DefaultMCMC::acceptance_probability ( const std::vector<double>& current_theta, const std::vector<double>& new_theta ) {
	double qnew;
	unsigned int i;
	qnew     = - getModel()->neglpost ( new_theta, getData() );
	for (i=0; i<getModel()->getNparams(); i++) {
		qnew -= log ( proposaldistributions[i]->pdf ( new_theta[i] ) );
	}

	/*
	std::cerr << "qnew = " << qnew << "\n"
		<< "Points: " << std::setiosflags(std::ios::fixed) << new_theta[0] << " " << new_theta[1] << " " << new_theta[2] << "\n"
		<< "Q*:     " << std::setiosflags(std::ios::fixed) << log ( proposaldistributions[0]->pdf ( new_theta[0] ) )
				<< " " << log ( proposaldistributions[1]->pdf ( new_theta[1] ) )
				<< " " << log ( proposaldistributions[2]->pdf ( new_theta[2] ) )
//				<< " " << log ( proposaldistributions[3]->pdf ( new_theta[3] ) )
				<< "\n"
		<< "P*:     " << -getModel()->neglpost ( new_theta, getData () ) << "\n"
		<< "qnew = " << qnew << "\n"
		<< "qold = " << qold << "\n"
		<< "qnew-qold = " << qnew - qold << " P(acceptance) = " << exp(qnew-qold) << "\n";
	*/

#ifdef DEBUG_MCMC
	std::cerr << "p(accept) = " << exp(qnew-qold) << "\n";
#endif

	return qnew;
}

void DefaultMCMC::proposePoint (
		std::vector<double> &current_theta,
		std::vector<double> &stepwidths,
		PsiRandom * proposal,
		std::vector<double> &new_theta ) {
	unsigned int i;

	for ( i=0; i<new_theta.size(); i++ ) {
		new_theta[i] = proposaldistributions[i]->rand ();
#ifdef DEBUG_MCMC
        std::cerr << new_theta[i] << "\t";
#endif
	}
#ifdef DEBUG_MCMC
	std::cerr << "\n";
#endif
}

/**********************************************************************
 *
 * Hybird MCMC
 *
 */

HybridMCMC::HybridMCMC ( const PsiPsychometric* Model, const PsiData* Data, int Nleap )
	: PsiSampler ( Model, Data ),
	currenttheta ( Model->getStart( Data ) ),
	newtheta     ( Model->getNparams(), 0 ),
	momentum  ( Model->getNparams(), 0 ),
	gradient (     Model->getNparams(), 0 ),
	currentgradient ( Model->getNparams(), 0 ),
	stepsizes (    Model->getNparams(), .001),
	Nleapfrog(Nleap),
	Naccepted (0)
{
	proposal = new GaussRandom;

	setTheta ( currenttheta );

	stepsizes[0] = 0.001;
	stepsizes[1] = 0.001;
	stepsizes[2] = 0.0001;
}

std::vector<double> HybridMCMC::draw ( void ) {
	unsigned int i;
	const PsiPsychometric * model ( getModel() );
	const PsiData *         data  ( getData() );

	for (i=0; i<model->getNparams(); i++)
		momentum[i] = proposal->draw();

	currentH = 0;
	for (i=0; i<model->getNparams(); i++)
		currentH += momentum[i] * momentum[i];
	currentH *= 0.5;
	currentH += energy;

	leapfrog();

	newenergy = model->neglpost ( newtheta, data );
	newH = 0;
	for (i=0; i<model->getNparams(); i++)
		newH += momentum[i] * momentum[i];
	newH *= 0.5;
	newH += newenergy;

	if ( log(proposal->rngcall()) < currentH-newH ) {
		// Accept
		for (i=0; i<model->getNparams(); i++) {
			currenttheta[i] = newtheta[i];
			currentgradient[i] = gradient[i];
		}
		energy = newenergy;
		Naccepted ++;
#ifdef DEBUG_MCMC
		std::cerr << " * ";
#endif
	}
#ifdef DEBUG_MCMC
	else
		std::cerr << "   ";
	std::cout << currenttheta[0] << "\n";
#endif
	return currenttheta;
}

void HybridMCMC::setTheta ( const std::vector<double>& theta ) {
	unsigned int i;
	currenttheta = theta;

	for (i=0; i<getModel()->getNparams(); i++) {
		gradient[i] = getModel()->dlposteri ( currenttheta, getData(), i );
	}
	energy = getModel()->neglpost ( currenttheta, getData() );
}

void HybridMCMC::setStepSize ( const std::vector<double>& sizes ) {
	if (sizes.size()==stepsizes.size())
		stepsizes = sizes;
	else
		throw BadArgumentError();
}

void HybridMCMC::setStepSize ( double size, unsigned int param ) {
	if ( param>=stepsizes.size() )
		throw BadIndexError();

	stepsizes[param] = size;
}

void HybridMCMC::leapfrog ( void ) {
	int i,n;
	int Nparams(getModel()->getNparams());
	const PsiPsychometric * model (getModel());

	gradient = currentgradient;
	newtheta = currenttheta;

	for (n=0; n<Nleapfrog; n++) {
		for (i=0; i<Nparams; i++)
			momentum[i] -= 0.5 * stepsizes[i] * gradient[i];

		for (i=0; i<Nparams; i++)
			newtheta[i] +=          stepsizes[i] * momentum[i];

		for (i=0; i<Nparams; i++)
			gradient[i] = model->dlposteri ( newtheta, getData(), i );

		for (i=0; i<Nparams; i++)
			momentum[i] -= 0.5 * stepsizes[i] * gradient[i];
	}
}

double HybridMCMC::getDeviance ( void ) {
	return getModel()->deviance ( currenttheta, getData() );
}

MCMCList HybridMCMC::sample ( unsigned int N ) {
	MCMCList out ( N, getModel()->getNparams(), getData()->getNblocks() );
	unsigned int i;

	for (i=0; i<N; i++) {
		out.setEst ( i, draw(), 0. );
		out.setdeviance ( i, getDeviance() );
#ifdef DEBUG_MCMC
		std::cerr
			<< "H: "       << std::setiosflags(std::ios::fixed) << currentH
			<< " H_: "     << std::setiosflags(std::ios::fixed) << newH
			<< " E: "      << std::setiosflags(std::ios::fixed) << energy
			<< " accept: " << std::setiosflags(std::ios::fixed) << double(Naccepted)/(i+1) << "\n";
#endif
	}

#ifdef DEBUG_MCMC
	std::cerr << "Acceptance rate: " << double(Naccepted)/N << "\n";
#endif

	return out;
}

/**********************************************************************
 *
 * Evidence
 *
 */

double ModelEvidence ( const PsiPsychometric* pmf, const PsiData* data )
{
	std::vector<double> prm ( pmf->getNparams() );
	double E(0);
	unsigned int i,k,n(50000);

	for ( i=0; i<n; i++ ) {
		for ( k=0; k<pmf->getNparams(); k++ )
			prm[k] = pmf->randPrior ( k );

		E += exp ( - pmf->negllikeli ( prm, data ) );
	}

	E /= n;

	return E;
}

std::vector<double> OutlierDetection ( const PsiPsychometric* pmf, OutlierModel* outl, const PsiData* data )
{
	unsigned int i;
	std::vector<double> out ( data->getNblocks() );
	double E ( ModelEvidence ( pmf, data ) );

	for ( i=0; i<data->getNblocks(); i++ ) {
		outl->setexclude ( i );
		out[i] = E/ModelEvidence ( outl, data );
	}

	return out;
}
