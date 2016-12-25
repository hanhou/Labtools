#include "getstart.h"

std::vector<double> linspace ( double xmin, double xmax, unsigned int n ) {
	double dummy;
	if ( xmin>xmax ) { // make sure that xmin is really less than xmax
		dummy = xmin;
		xmin = xmax;
		xmax = dummy;
	}
	unsigned int i;
	double xstep ( (xmax-xmin)/(n-1) );
	std::vector<double> out (n);

	out[0] = xmin;
	for (i=1; i<n; i++) {
		out[i] = out[i-1] + xstep;
	}

	return out;
}

/************************************************** PsiGrid methods ********************************************/

PsiGrid::PsiGrid ( const std::vector<double>& xmin, const std::vector<double>& xmax, unsigned int gridsize ) :
	ndim ( xmin.size() ),
	ngrid ( gridsize ),
	grid1d ( xmin.size() ),
	lower_bounds(xmin),
	upper_bounds(xmax)
{
	if ( lower_bounds.size() != upper_bounds.size() )
		throw PsiError ( "Upper and lower grid bounds are unequal" );
	unsigned int i;

	for ( i=0; i<ndim; i++ ) {
		grid1d[i] = linspace ( lower_bounds[i], upper_bounds[i], gridsize );
	}
}

PsiGrid PsiGrid::shift ( const std::vector<double>& newposition ) const
{
	std::vector<double> xmin ( lower_bounds );
	std::vector<double> xmax ( upper_bounds );
	unsigned int i;
	double gridmid;

	for ( i=0; i<newposition.size (); i++ ) {
		gridmid = (xmax[i]-xmin[i])/2.;
		xmin[i] += newposition[i]-gridmid;
		xmax[i] += newposition[i]-gridmid;
	}

	return PsiGrid( xmin, xmax, get_gridsize() );
}

PsiGrid PsiGrid::shrink ( const std::vector<double>& newposition ) const
{
	std::vector<double> xmin ( lower_bounds );
	std::vector<double> xmax ( upper_bounds );
	unsigned int i;
	double xstep;

	for ( i=0; i<newposition.size(); i++ )
	{
		xstep = grid1d[i][1]-grid1d[i][0];
		xmin[i] = newposition[i]-xstep;
		xmax[i] = newposition[i]+xstep;
	}

	return PsiGrid ( xmin, xmax, get_gridsize() );
}

PsiGrid PsiGrid::subgrid ( void ) const
{
	std::vector<double> xmin ( lower_bounds.size()-1 );
	std::vector<double> xmax ( upper_bounds.size()-1 );
	unsigned int i;

	for ( i=0; i<xmin.size(); i++ ) {
		xmin[i] = lower_bounds[i+1];
		xmax[i] = upper_bounds[i+1];
	}

	return PsiGrid ( xmin, xmax, get_gridsize() );
}

/************************************************** PsiGrid functions ********************************************/

void makegridpoints (
		const PsiGrid& grid,
		std::vector<double> prm,
		unsigned int pos,
		std::list< std::vector<double> > *gridpoints
		)
{
	unsigned int i;
	if ( grid.dimension() != prm.size() ) {
		throw PsiError ( "grid and parameter vector don't match" );
	}

	if ( pos>=grid.dimension() ) {
		// We are in the lowest level of iteration
		gridpoints->push_back ( prm );
		return;
	} else {
		// We have to loop over this level
		for ( i=0; i<grid.get_gridsize(); i++ ) {
			prm[pos] = grid(pos,i);
			makegridpoints ( grid, prm, pos+1, gridpoints );
		}
	}
}

std::vector<double> pymakegridpoints (
		const PsiGrid& grid,
		std::vector<double> prm,
		unsigned int pos
		)
{
	std::list < std::vector<double> >gridpoints;
	std::list < std::vector<double> >::const_iterator griditer;
	makegridpoints ( grid, prm, pos, &gridpoints );
	griditer = gridpoints.begin();
	unsigned int nparams = griditer->size();
	unsigned int npoints = gridpoints.size();
	std::cerr << "Gridpoints:" << npoints << "\nParams:" << nparams << "\n";

	std::vector<double> out ( nparams*npoints );
	unsigned int i,j;

	for (griditer=gridpoints.begin(),i=0; griditer!=gridpoints.end(); i+= nparams, griditer++) {
		for (j=0; j<nparams; j++)
			out[i+j] = griditer->at(j);
	}
	return out;
}


void evalgridpoints (
		const std::list< std::vector<double> >& gridpoints,
		std::list< std::vector<double> > *bestprm,
		std::list< double > *L,
		const PsiData* data,
		const PsiPsychometric* pmf,
		unsigned int nbest
		)
{
	std::list< std::vector<double> >::const_iterator griditer;
	std::list< std::vector<double> >::iterator iter_prm;
	std::list< double >::iterator iter_L;
	double l;
	double a,b;
	std::vector<double> prm;
	const PsiCore *core = pmf->getCore();
	bool store(true);

	for ( griditer=gridpoints.begin(); griditer!=gridpoints.end(); griditer++ ) {
		// Transform parameters and get negative log posterior
		a = (*griditer)[0];
		b = (*griditer)[1];
		b = 1./b;
		a = -a*b;
		// prm = core->transform ( pmf->getNparams(), 1./b, -a/b );
		prm = core->transform ( pmf->getNparams(), a, b );
		prm[2] = (*griditer)[2];
		if ( pmf->getNparams() > 3 ) prm[3] = (*griditer)[3];
		l = pmf->neglpost ( prm, data );

		// Where does it belong?
		for ( iter_L=L->begin(), iter_prm=bestprm->begin() ; iter_L!=L->end(); iter_L++, iter_prm++ ) {
			if ( l==(*iter_L) ) {
				if ( (*iter_prm) == (*griditer) )
					store = false;
				else
					store = true;
				break;
			} else if ( l<(*iter_L) ) {
				store = true;
				break;
			} else {
				store = false;
			}
		}
		// insert the values if they are good enough
		if ( store ) {
			L->insert ( iter_L, l );
			bestprm->insert ( iter_prm, std::vector<double>(*griditer) );
		}

		// Reduce the size of the best parameters list
		while ( L->size() > nbest ) {
			L->pop_back();
			bestprm->pop_back();
		}
	}
}

void updategridpoints (
		const PsiGrid& grid,
		const std::list< std::vector<double> >& bestprm,
		std::list< std::vector<double> > *newgridpoints,
		std::list< PsiGrid > *newgrids
		)
{
	// modify grid size: If a point in bestprm is on the edge of the grid: make a new, larger grid
	// if a point is an interior point of the grid, shring the grid to the area around that grid

	std::list< std::vector<double> >::const_iterator iter_prm;
	std::vector<double> prm ( bestprm.front().size() );
	bool isedge (false);
	unsigned int i;
	PsiGrid newgrid;

	for ( iter_prm=bestprm.begin(); iter_prm!=bestprm.end(); iter_prm++ ) {
		// Check whether the current point is on the edge of the grid
		isedge = false;
		for ( i=0; i<iter_prm->size(); i++ ) {
			isedge += (*iter_prm)[i]==grid.get_lower(i);
			isedge += (*iter_prm)[i]==grid.get_upper(i);
		}

		if (isedge) {
			newgrid = grid.shift ( *iter_prm );
		} else {
			newgrid = grid.shrink ( *iter_prm);
		}
		makegridpoints ( newgrid, prm, 0, newgridpoints );
		newgrids->push_back ( newgrid );
	}
}

/*************************************** Range heuristics ****************************************/

void a_range ( const PsiData* data, double *xmin, double *xmax ) {
	double x;
	unsigned int i;
	*xmin = 1e5;
	*xmax = -1e5;

	// Heuristic:
	// a will be between lowest and highest stimulus level
	for ( i=0; i<data->getNblocks(); i++ ) {
		x = data->getIntensity ( i );
		if ( x<*xmin ) {
			*xmin = x;
		}
		if ( x>*xmax ) {
			*xmax = x;
		}
	}
}

void b_range ( const PsiData* data, double *xmin, double *xmax ) {
	double x,p(1),xx,pp(0),pc;
	std::vector<double> intensities ( data->getIntensities() );
	unsigned int i,j;
	double d;
	*xmin = 1e5;
	*xmax = -1e5;

	// Heuristic:
	// b will be between rising over the whole range of stimulus levels and
	// rising between two adjacent stimulus levels

	// First determine step sizes
	for ( i=0; i<intensities.size(); i++ ) {
		for ( j=i; j<intensities.size(); j++ ) {
			d = fabs ( intensities[i] - intensities[j] );
			if (d==0) continue;

			if (d>*xmax) {
				*xmax = d;
			}

			if (d<*xmin) {
				*xmin = d;
			}
		}
	}

	// Is the psychometric function rising or falling overall
	for ( i=0; i<intensities.size(); i++ ) {
		pc = data->getPcorrect ( i );
		if ( pc<p ) {
			p = pc;
			x = intensities[i];
		}
		if ( pc>pp ) {
			pp = pc;
			xx = intensities[i];
		}
	}
	if ( xx<x ) { // psychometric function is falling
		x = *xmin;
		*xmin = *xmax;
		*xmax = x;
	}

	// In any case, b should be scaled to be roughly equivalent to w
	/*
	x = 2*log(9.);
	*xmin /= x;
	*xmax /= x;
	*/
}

void lm_range ( const PsiData* data, double *xmin, double *xmax ) {
	double p,pmax(0);
	unsigned int i;

	// Heuristic:
	// lm goes from 0 to twice the distance between 1 and the highest response probability
	*xmin = 0;
	for ( i=0; i<data->getNblocks(); i++ ) {
		p = data->getPcorrect ( i );
		if ( p > pmax ) {
			pmax = p;
		}
	}
	*xmax = 2*(1-pmax);
	if (*xmax>=1) *xmax=.99;
	if (*xmax<=0) *xmax=0.1;
}

void gm_range ( const PsiData* data, double *xmin, double *xmax ) {
	double p, pmin(0);
	unsigned int i;

	// Heuristic:
	// gm goes from 0 to twice the lowsest response probability
	*xmin = 0;
	for ( i=0; i<data->getNblocks(); i++ ) {
		p = data->getPcorrect ( i );
		if ( p<pmin ) {
			pmin = p;
		}
	}

	*xmax = 2*pmin;
	if (*xmax>1) *xmax=.99;
	if (*xmax<0.1) *xmax=.1;
}

void parameter_range ( const PsiData* data, const PsiPsychometric * pmf, unsigned int prmindex, double *xmin, double *xmax ) {
	// Call the correct initial range function for the parameter
	double mps,mms,s;
	const PsiPrior *prior = pmf->getPrior( prmindex );
	// incorporate information from prior
	mps = mms = prior->mean();
	s = prior->std();
	// mean plus (+) std-deviation
	mps += 2*s;
	// mean minus (-) std-deviation
	mms -= 2*s;

	switch ( prmindex ) {
		case 0:
			a_range ( data, xmin, xmax );
			break;
		case 1:
			b_range ( data, xmin, xmax );
			break;
		case 2:
			lm_range ( data, xmin, xmax );
			break;
		case 3:
			gm_range ( data, xmin, xmax );
			break;
	}

	if ( *xmin < mms ) *xmin = mms;
	if ( *xmax > mps ) *xmax = mps;
}

std::vector<double> getstart (
		const PsiPsychometric* pmf,
		const PsiData* data,
		unsigned int gridsize,
		unsigned int nneighborhoods,
		unsigned int niterations,
		std::vector<double> *incr )
{
	std::vector<double> xmin ( pmf->getNparams() );
	std::vector<double> xmax ( pmf->getNparams() );
	std::list< std::vector<double> > bestprm;
	std::list< double > L;
	unsigned int i,j, ngrids;

	// Set up the initial grid
	for ( i=0; i<pmf->getNparams(); i++ ) {
		parameter_range ( data, pmf, i, &(xmin[i]), &(xmax[i]) );
	}

	// Make an initial grid
	PsiGrid grid ( xmin, xmax, gridsize ), currentgrid;
	std::list< PsiGrid > newgrids;
	newgrids.push_back ( grid );

	// Perform first evaluation on the grid
	std::list< std::vector<double> > gridpoints;
	makegridpoints ( grid, xmin, 0, &gridpoints );
	evalgridpoints ( gridpoints, &bestprm, &L, data, pmf, nneighborhoods );

	// potentially more evaluations
	for ( i=0; i<niterations; i++ ) {
		while ( newgrids.size() > nneighborhoods ) {
			newgrids.pop_front ();
		}

		ngrids = newgrids.size();
		for ( j=0; j<ngrids; j++ ) {
			currentgrid = newgrids.front ();
			newgrids.pop_front ();
			gridpoints = std::list< std::vector<double> > ();
			updategridpoints ( currentgrid, bestprm, &gridpoints, &newgrids );
			evalgridpoints ( gridpoints, &bestprm, &L, data, pmf, nneighborhoods );
		}

		// evalgridpoints ( gridpoints, &bestprm, &L, data, pmf, nneighborhoods );
	}

	// Now transform the best parameter to the suitable format
	const PsiCore *core = pmf->getCore();
	double a ( bestprm.front()[0] ), b ( bestprm.front()[1] );
	// std::cerr << "Raw starting values:";
	// for (i=0; i<bestprm.front().size(); i++) {
	// 	std::cerr << " " << bestprm.front()[i];
	// }
	// std::cerr << "\n";
	b = 1./b;
	a = -a*b;
	std::vector<double> out = core->transform ( pmf->getNparams(), a, b );
	out[2] = bestprm.front()[2];
	if ( pmf->getNparams() > 3 ) out[3] = bestprm.front()[3];

	if ( incr!=NULL ) {
		if ( incr->size() != pmf->getNparams() ) throw ( BadArgumentError ( "Wrong size for incr" ) );
		currentgrid = newgrids.front();
		for ( i=0; i<pmf->getNparams(); i++ ) {
			(*incr)[i] = 10*currentgrid.get_incr(i);
		}
	}

	// std::cerr << "Starting values:";
	// for (i=0; i<out.size(); i++)
	// 	std::cerr << " " << out[i];
	// std::cerr << "\n";

	return out;
}
