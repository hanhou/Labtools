/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#include "data.h"

/************************************************************
 * Constructors                                             *
 ************************************************************/
PsiData::PsiData (
	std::vector<double> x,
	std::vector<int>    N,
	std::vector<int>    k,
	int nAFC
	) :
	intensities(x), Ntrials(N), Ncorrect(k), Pcorrect(k.size()), logNoverK(k.size()), Nalternatives(nAFC)
{
	unsigned int i,n;
	for ( i=0; i<k.size(); i++ ) {
		Pcorrect[i] = double(Ncorrect[i])/Ntrials[i];
		logNoverK[i] = 0;
		for ( n=1; n<=(unsigned int) (k[i]); n++ )
			logNoverK[i] += log(N[i]+1-n) - log(n);
	}
}

PsiData::PsiData (
	std::vector<double> x,
	std::vector<int>    N,
	std::vector<double> p,
	int nAFC
	) :
	intensities(x), Ntrials(N), Ncorrect(p.size()), Pcorrect(p), Nalternatives(nAFC)
{
	unsigned int i;
	double k;
	for ( i=0; i<p.size(); i++ ) {
		k = Ntrials[i] * Pcorrect[i];
		if ( fabs(k-int(k)) > 1e-7 )    // The fraction of correct responses does not correspond to an integer number of correct responses
			std::cerr << "WARNING: fraction of correct responses does not correspond to an integer number of correct responses!\n";
		Ncorrect[i] = int(k);
	}
}

/************************************************************
 * Setters                                                  *
 ************************************************************/

void PsiData::setNcorrect ( const std::vector<int>& newNcorrect )
{
	Ncorrect = newNcorrect;
	unsigned int i;
	for ( i=0; i<Ncorrect.size(); i++ )
		Pcorrect[i] = double(Ncorrect[i])/Ntrials[i];
}

/************************************************************
 * Getters                                                  *
 ************************************************************/

// TODO: implement these as inline functions?

const std::vector<double>& PsiData::getIntensities ( void ) const
{
    return intensities;
}

const std::vector<int>& PsiData::getNtrials ( void ) const
{
    return Ntrials;
}

const std::vector<int>& PsiData::getNcorrect ( void ) const
{
    return Ncorrect;
}

const std::vector<double>& PsiData::getPcorrect ( void ) const
{
    return Pcorrect;
}

double PsiData::getIntensity ( unsigned int i ) const
{
	if ( i>=0 && i<intensities.size() )
		return intensities[i];
	else
		throw BadIndexError();
}

int PsiData::getNtrials ( unsigned int i ) const
{
	if ( i>=0 && i<Ntrials.size() )
		return Ntrials[i];
	else
		throw BadIndexError();
}

int PsiData::getNcorrect ( unsigned int i ) const
{
	if ( i>=0 && i<Ncorrect.size() )
		return Ncorrect[i];
	else
		throw BadIndexError();
}

double PsiData::getPcorrect ( unsigned int i ) const
{
	if ( i>=0 && i<Pcorrect.size() )
		return Pcorrect[i];
	else
		throw BadIndexError();
}

double PsiData::getNoverK ( unsigned int i ) const
{
	if ( i>=0 && i<logNoverK.size() )
		return logNoverK[i];
	else
		throw BadIndexError();
}

int PsiData::getNalternatives ( void ) const
{
	return Nalternatives;
}

std::vector<int> PsiData::nonasymptotic ( void ) const
{
	unsigned int i,j,Ngood(0);
	double guess ( 1./Nalternatives );
	if ( Nalternatives<2 )
		guess = 0;

	for (i=0; i<getNblocks(); i++) {
		// if ( Pcorrect[i] < 0.95 && Pcorrect[i] > guess+0.05  )
		if ( Pcorrect[i] < 1. )
			Ngood ++;
	}

	std::vector<int> out ( Ngood );
	j=0;
	for (i=0; i<getNblocks(); i++) {
		// if ( Pcorrect[i] < 0.95 && Pcorrect[i] > guess+0.05 )
		if ( Pcorrect[i] < 1. )
			out[j++] = i;
	}

	return out;
}
