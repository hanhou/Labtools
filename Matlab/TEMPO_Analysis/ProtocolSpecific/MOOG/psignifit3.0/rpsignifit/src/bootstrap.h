/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H

#include <vector>
#include <cmath>
#include "psychometric.h"
#include "mclist.h"
#include "optimizer.h"

/** \brief perform a parametric bootstrap
 *
 * A parametric bootstrap is performed by sampling from a binomial distribution with success probability given by the psychometric
 * function. if BCa is true, bias correction and acceleration constant are calculated for the cuts given in cuts.
 */
BootstrapList bootstrap (
		unsigned int B,                        ///< number of bootstrap samples
		const PsiData * data,         ///< data that are to form the basis of the whole procedure
		const PsiPsychometric* model, ///< model to be fitted
		std::vector<double> cuts,     ///< performance levels at which the threshold should be calculated
		std::vector<double>* param=NULL,   ///< parameter vector on which parametric bootstrap should be based
		bool BCa=true,                ///< calculate bias correction and acceleration?
		bool parametric=true          ///< Perform parametric bootstrap?
		);

/** \brief perform jackkifing to detect influential observations and outliers
 *
 * Wichmann & Hill (2001) suggest performing jackknife resampling on the data to determine
 * influential observations and outliers.
 *
 * Wichmann & Hill (2001) The psychometric function: I. Fitting, sampling, and goodness of fit. Perception & Psychophysics, 63(8), 1293--1313.
 */
JackKnifeList jackknifedata ( const PsiData * data, const PsiPsychometric* model );


void newsample ( const PsiData * data, const std::vector<double>& p, std::vector<int> * sample );

#endif
