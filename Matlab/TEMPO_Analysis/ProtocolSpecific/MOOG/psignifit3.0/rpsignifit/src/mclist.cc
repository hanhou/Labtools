/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#include "mclist.h"

void newsample ( const PsiData * data, const std::vector<double>& p, std::vector<int> * sample ) {
	/* Draw a new sample from the psychometric function */
	BinomialRandom binomial ( 10, 0.5 );    // Initialize with nonsense parameters
	unsigned int k;                                            // Block index

	for ( k=0; k<data->getNblocks(); k++ ) {
		binomial.setprm ( data->getNtrials(k), p[k] );
		(*sample)[k] = binomial.draw ();
	}
}

/************************************************************
 * PsiMClist methods
 */

std::vector<double> PsiMClist::getEst ( unsigned int i ) const
{
	// Check that the call does not ask for something we don't have
	if ( i>=getNsamples() )
		throw BadIndexError();

	unsigned int k;
	std::vector<double> out ( getNparams() );

	for (k=0; k<getNparams(); k++)
		out[k] = mcestimates[k][i];

	return out;
}

double PsiMClist::getEst ( unsigned int i, unsigned int prm ) const
{
	// check that the call does not ask for something we don't have
	if ( i>=getNsamples() )
		throw BadIndexError();
	if ( prm>=getNparams() )
		throw BadIndexError();

	return mcestimates[prm][i];
}

void PsiMClist::setEst ( unsigned int i, const std::vector<double> est, double deviance )
{
	// Check that the call does not ask for something we can't do
	if ( i>=getNsamples() )
		throw BadIndexError();

	unsigned int k;
	for ( k=0; k<getNparams(); k++ )
		mcestimates[k][i] = est[k];
	deviances[i] = deviance;
}

double PsiMClist::getPercentile ( double p, unsigned int prm ) {
	if ( prm>=getNparams() )
		throw BadIndexError();
	if ( p>1 || p<0 )
		throw BadArgumentError();

	int position;
	sort( mcestimates[prm].begin(),mcestimates[prm].end() );

	position = getNsamples()*p;

	return mcestimates[prm][position];
}

void PsiMClist::setdeviance ( unsigned int i, double deviance ) {
	if ( i>=getNsamples() )
		throw BadIndexError();

	deviances[i] = deviance;
}

double PsiMClist::getdeviance ( unsigned int i ) const {
	if ( i>=getNsamples() )
		throw BadIndexError();

	return deviances[i];
}

double PsiMClist::getDeviancePercentile ( double p ) {
	if ( p<=0 || p>= 1 )
		throw BadArgumentError();

	int ind ( p*deviances.size() );

	sort( deviances.begin(), deviances.end() );

	return deviances[ind];
}

double PsiMClist::getMean ( unsigned int prm ) const {
	double m(0);
	unsigned int i,Nsamples(getNsamples());
	if ( prm>=getNparams() )
		throw BadIndexError();

	for (i=0; i<Nsamples; i++)
		m += getEst ( i, prm );

	m /= Nsamples;
	return m;
}

double PsiMClist::getStd ( unsigned int prm ) const {
	double m ( getMean ( prm ) ), s(0), ss;
	unsigned int i, Nsamples(getNsamples());
	if ( prm>=getNparams() )
		throw BadIndexError();

	for (i=0; i<Nsamples; i++) {
		ss = getEst ( i, prm ) - m;
		s += ss*ss;
	}

	s /= Nsamples-1;
	return sqrt(s);
}

/************************************************************
 * BootstrapList methods
 */

void BootstrapList::setData ( unsigned int i, const std::vector<int>& newdata )
{
	if ( i>=getNsamples() || i<0 )
		throw BadIndexError();

	unsigned int k;
	for ( k=0; k<getNblocks(); k++ )
		data[i][k] = newdata[k];
}

std::vector<int> BootstrapList::getData ( unsigned int i ) const
{
	if ( i>=getNsamples() || i<0 )
		throw BadIndexError();

	return data[i];
}

double BootstrapList::getThres ( double p, unsigned int cut ) {
	if ( cut>=cuts.size() )
		throw BadIndexError();
	if ( p>1 || p<0 )
		throw BadArgumentError();

	int position;
	sort( thresholds[cut].begin(), thresholds[cut].end() );

	// Bias correction of p
	if (BCa)
		p = Phi(bias_t[cut] + (invPhi(p) + bias_t[cut])/(1-acceleration_t[cut]*(invPhi(p) + bias_t[cut])));

	position = int(getNsamples()*p);

	return thresholds[cut][position];
}

double BootstrapList::getThres_byPos ( unsigned int i, unsigned int cut ) {
	if ( cut>=cuts.size() )
		throw BadIndexError();
	if (i>getNsamples())
		throw BadIndexError();

	return thresholds[cut][i];
}

void BootstrapList::setThres ( double thres, unsigned int i, unsigned int cut )
{
	if ( i>=getNsamples() )
		throw BadIndexError();
	if (cut>=cuts.size() )
		throw BadIndexError();

	thresholds[cut][i] = thres;
}

double BootstrapList::getSlope ( double p, unsigned int cut ) {
	if ( cut>=cuts.size() )
		throw BadIndexError();
	if ( p>1 || p<0 )
		throw BadArgumentError();

	int position;
	sort( slopes[cut].begin(), slopes[cut].end() );

	// Bias correction of p
	if (BCa)
		p = Phi(bias_s[cut] + (invPhi(p) + bias_s[cut])/(1-acceleration_s[cut]*(invPhi(p) + bias_s[cut])));

	position = int(getNsamples()*p);

	return slopes[cut][position];
}

double BootstrapList::getSlope_byPos ( unsigned int i, unsigned int cut ) {
	if ( cut>=cuts.size() )
		throw BadIndexError();
	if (i>getNsamples())
		throw BadIndexError();

	return slopes[cut][i];
}

void BootstrapList::setSlope ( double slope, unsigned int i, unsigned int cut )
{
	if ( i>=getNsamples() )
		throw BadIndexError();
	if (cut>=cuts.size() )
		throw BadIndexError();

	slopes[cut][i] = slope;
}

double BootstrapList::getCut ( unsigned int i ) const
{
	if ( i>=cuts.size() || i<0 )
		throw BadIndexError();

	return cuts[i];
}

void BootstrapList::setRpd ( unsigned int i, double r_pd ) {
	if ( i>=getNsamples() )
		throw BadIndexError();

	Rpd[i] = r_pd;
}

double BootstrapList::getRpd ( unsigned int i ) const {
	if ( i>=getNsamples() )
		throw BadIndexError();

	return Rpd[i];
}

double BootstrapList::percRpd ( double p ) {
	if ( p<0 || p>1 )
		throw BadArgumentError();

	int index ( p*(getNsamples()-1));

	sort ( Rpd.begin(), Rpd.end() );

	return Rpd[index];
}

void BootstrapList::setRkd ( unsigned int i, double r_kd ) {
	if ( i>=getNsamples() )
		throw BadIndexError();

	Rkd[i] = r_kd;
}

double BootstrapList::getRkd ( unsigned int i ) const {
	if ( i>=getNsamples() )
		throw BadIndexError();

	return Rkd[i];
}

double BootstrapList::percRkd ( double p ) {
	if ( p<0 || p>1 )
		throw BadIndexError();

	int index ( p*(getNsamples()-1) );

	sort ( Rkd.begin(), Rkd.end() );

	return Rkd[index];
}

/************************************************************
 * JackKnifeList methods
 */

double JackKnifeList::influential ( unsigned int block, const std::vector<double>& ci_lower, const std::vector<double>& ci_upper ) const {
	unsigned int prm;
	double est;
	double infl(0),x;

	for ( prm=0; prm<getNparams(); prm++ ) {
		est = getEst(block,prm);
		if ( est<mlestimate[prm] ) {
			x = ( mlestimate[prm]-est ) / ( mlestimate[prm]-ci_lower[prm] );
		} else {
			x = ( est-mlestimate[prm] ) / ( ci_upper[prm]-mlestimate[prm] );
		}
		if ( x>infl )
			infl = x;
	}
	return infl;
}

bool JackKnifeList::outlier ( unsigned int block ) const {
	if ( block>=getNblocks() )
		throw BadIndexError();

	if ( maxdeviance-getdeviance(block) > 6.63 )
		return true;
	return false;
}

/************************************************************
 * MCMCList methods
 */

void MCMCList::setppData ( unsigned int i, const std::vector<int>& ppdata, double ppdeviance )
{
	if ( i>=getNsamples() || i<0 )
		throw BadIndexError ();

	unsigned int k;
	for ( k=0; k<getNblocks(); k++ )
		posterior_predictive_data[i][k] = ppdata[k];
	posterior_predictive_deviances[i] = ppdeviance;
}

std::vector<int> MCMCList::getppData ( unsigned int i ) const
{
	if ( i>=getNsamples() || i<0 )
		throw BadIndexError();

	return posterior_predictive_data[i];
}

int MCMCList::getppData ( unsigned int i, unsigned int j ) const
{
	if ( i>=getNsamples() )
		throw BadIndexError();
	if ( j>=getNblocks() )
		throw BadIndexError();

	return posterior_predictive_data[i][j];
}

double MCMCList::getppDeviance ( unsigned int i ) const
{
	if ( i>=getNsamples() )
		throw BadIndexError();

	return posterior_predictive_deviances[i];
}

void MCMCList::setppRpd ( unsigned int i, double Rpd )
{
	if ( i>=getNsamples() )
		throw BadIndexError();

	posterior_predictive_Rpd[i] = Rpd;
}

double MCMCList::getppRpd ( unsigned int i ) const
{
	if ( i>=getNsamples() )
		throw BadIndexError();

	return posterior_predictive_Rpd[i];
}

void MCMCList::setppRkd ( unsigned int i, double Rkd )
{
	if ( i>=getNsamples() )
		throw BadIndexError();

	posterior_predictive_Rkd[i] = Rkd;
}

double MCMCList::getppRkd ( unsigned int i ) const
{
	if ( i>=getNsamples() )
		throw BadIndexError();

	return posterior_predictive_Rkd[i];
}

void MCMCList::setRpd ( unsigned int i, double Rpd )
{
	if ( i>=getNsamples() )
		throw BadIndexError();

	posterior_Rpd[i] = Rpd;
}

double MCMCList::getRpd ( unsigned int i ) const
{
	if ( i>=getNsamples() )
		throw BadIndexError();

	return posterior_Rpd[i];
}

void MCMCList::setRkd ( unsigned int i, double Rkd )
{
	if ( i>=getNsamples() )
		throw BadIndexError();

	posterior_Rkd[i] = Rkd;
}

double MCMCList::getRkd ( unsigned int i ) const
{
	if ( i>=getNsamples() )
		throw BadIndexError();

	return posterior_Rkd[i];
}


void MCMCList::setlogratio ( unsigned int i, unsigned int j, double logratio )
{
	if ( i>=getNsamples() )
		throw BadIndexError();
	if ( j>=getNblocks() )
		throw BadIndexError();

	logratios[i][j] = logratio;
}

double MCMCList::getlogratio ( unsigned int i, unsigned int j ) const
{
	if ( i>=getNsamples() )
		throw BadIndexError();
	if ( j>=getNblocks() )
		throw BadIndexError();

	return logratios[i][j];
}
