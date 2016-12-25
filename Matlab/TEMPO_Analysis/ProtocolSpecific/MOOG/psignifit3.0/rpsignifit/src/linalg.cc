/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#include "linalg.h"
#include "limits.h"

double sign ( double x ) {
	return x/fabs(x);
}

double househ ( const std::vector<double> *x, std::vector<double> *u ) {
	unsigned int i;
	double h;

	h = 0;
	for ( i=0; i<x->size(); i++ ) {
		h += (*x)[i]*(*x)[i];
		(*u)[i] = (*x)[i];
	}
	h = sqrt(h);

	if ( (*x)[0] == 0 ) (*u)[0] = h;
	else                (*u)[0] = (*x)[0] + sign((*x)[0])*h;

	for ( i=u->size()-1; i<UINT_MAX; i-- ) { // This loop should stop after (!) i has reached 0. As i is unsigned this will mean that i is MAXINT
		(*u)[i] /= (*u)[0];
	}

	return 1+fabs((*x)[0])/h;
}

double uuA ( const std::vector<double> *u, const Matrix *A, int j, int k, int l ) {
	unsigned int i, ii;
	std::vector<double> uA ( A->getncols() - j, 0 );
	for ( i=0; i<uA.size(); i++ )
		for ( ii=0; ii<u->size(); ii++ )
			uA[i] += (*u)[ii] * (*A)(j+ii,j+i);
	return (*u)[k] * uA[l];
}

Matrix::Matrix ( const std::vector< std::vector<double> >& A )
	: nrows(A.size()), ncols(A[0].size())
{
	data = new double [nrows*ncols];

	unsigned int i,j;
	for (i=0; i<nrows; i++) {
		for (j=0; j<ncols; j++) {
			this->operator() (i,j) = A[i][j];
		}
	}
}

Matrix::Matrix ( unsigned int nrows, unsigned int ncols)
	: nrows(nrows), ncols(ncols)
{
	data = new double [nrows*ncols];
	unsigned int i;
	for (i=0; i<nrows*ncols; i++)
		data[i] = 0;
}

Matrix::Matrix ( const Matrix& A )
	: nrows(A.getnrows()), ncols(A.getncols())
{
	data = new double [nrows*ncols];

	unsigned int i,j;
	for (i=0; i<nrows; i++) {
		for (j=0; j<ncols; j++) {
			this->operator() (i,j) = A(i,j);
		}
	}
}

double& Matrix::operator() ( unsigned int row, unsigned int col ) const
{
	//row --; col --;
	if (row>=nrows || col>=ncols) {
		throw MatrixError();
	}
	return data[row+nrows*col];
}

void Matrix::print ( void ) {
	unsigned int i,j;
	std::cout << "[ ";
	for (i=0; i<nrows; i++) {
		std::cout << "[";
		for (j=0; j<ncols; j++) {
			std::cout << " " << std::setprecision(3) << std::setw(5) << data[i+nrows*j] << (j!=ncols-1 ? " , " : (i==nrows-1 ? "] ]\n" : "],\n  ") );
		}
	}
}

Matrix* Matrix::cholesky_dec ( void ) const {
	if (nrows!=ncols)
		throw MatrixError();

	Matrix *L = new Matrix (nrows,nrows);
	unsigned int k,j,K;

	for (K=0; K<nrows; K++) {
		// l_{KK} = (a_{KK} - \sum_{k=1}^{K-1} l_{KK}^2)^{1/2}
		(*L)(K,K) = this->operator() (K,K);
		for (k=0; k<K; k++)
			(*L)(K,K) -= (*L)(K,k) * (*L)(K,k);
		(*L)(K,K) = sqrt((*L)(K,K));

		// l_{jK} = (a_{jK} - \sum_{k=1}^{K-1} l_{jk} l_{Kk})/l_{KK}, j >= K+1
		for (j=K+1; j<nrows; j++) {
			(*L)(j,K) = this->operator() (j,K);
			for (k=0; k<K; k++) {
				(*L)(j,K) -= (*L)(j,k)* (*L)(K,k);
			}
			(*L)(j,K) /= (*L)(K,K);
		}
	}

	return L;
}

Matrix* Matrix::lu_dec ( void ) const {
	if (nrows!=ncols) {
		throw MatrixError ();
	}

	Matrix *LU = new Matrix (*this);

	unsigned int i,j,k;
	double c;
	unsigned int pivotindex;
	double pivot;

	for (i=0; i<nrows-1; i++) {
		// Search Pivot element
		pivot = (*LU)(i,i);
		pivotindex = i;
		for (k=i+1; k<nrows; k++) {
			if ( fabs((*LU)(k,i))>pivot ) {
				pivot = fabs((*LU)(k,i));
				pivotindex = k;
			}
		}
		// Check that the pivot element does not vanish
		if ( pivot<1e-8 ) {
			delete LU;
			throw std::string ( "Matrix is numerically singular" );
		}
		// Swap pivot elements
		for (j=i; j<ncols; j++) {
			pivot = (*LU)(pivotindex,j);
			(*LU)(pivotindex,j) = (*LU)(i,j);
			(*LU)(i,j) = pivot;
		}

		// Eliminate
		for (k=i+1; k<nrows; k++) {
			c = (*LU)(k,i) / (*LU)(i,i);
			(*LU)(k,i) = c;
			for (j=i+1; j<nrows; j++) {
				(*LU)(k,j) = (*LU)(k,j) - c * (*LU)(i,j);
			}
		}
	}


	return LU;
}

Matrix* Matrix::qr_dec ( void ) const {
	Matrix *A = new Matrix ( *this );

	// A->print();

	int i, j, k, l;
	int m ( A->getnrows() ), n ( A->getncols() );
	std::vector<double> *u, *x;
	Matrix *uuAmat;
	double c;
	int min ( m-1>n ? n : m-1 );

	for ( j=0; j<min; j++ ) {

		x = new std::vector<double> ( m-j );
		u = new std::vector<double> ( m-j );
		uuAmat = new Matrix ( m-j,n-j );

		for ( i=j; i<m; i++ ) {
			(*x)[i-j] = (*A)(i,j);
		}
		c = househ ( x, u );

		for ( k=j; k<m; k++ ) {
			for ( l=j; l<n; l++ ) {
				(*uuAmat)(k-j,l-j) = uuA(u,A,j,k-j,l-j);
			}
	}

		for ( k=j; k<m; k++ ) {
			for ( l=j; l<n; l++ ) {
				(*A)(k,l) -= c * (*uuAmat)(k-j,l-j);
			}
		}

		delete x;
		delete u;
		delete uuAmat;
	}

	return A;
}

Matrix * Matrix::regularized_inverse ( double alpha ) const {
	// Tichonov regularized inverse
	if (nrows!=ncols)
		throw MatrixError ();
	int N(getnrows());

	Matrix *AAT = new Matrix (N,N);
	Matrix *rInv = new Matrix (N,N);
	std::vector<double> ATb (N);
	std::vector<double> x (N);

	int i,j,k;

	// Compute AAT
	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
			(*AAT)(i,j) = 0;
			for ( k=0; k<N; k++ )
				(*AAT)(i,j) += (*this)(i,k)*(*this)(k,j);
		}
	}

	// Regularize
	for (k=0; k<N; k++) (*AAT)(k,k) += alpha;

	// Solve for unit vectors
	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
			ATb[j] = (*this)(j,i);
		}
		x = AAT->solve(ATb);
		for (j=0; j<N; j++) {
			(*rInv)(i,j) = x[j];
		}
	}

	delete AAT;

	return rInv;
}

Matrix *Matrix::inverse_qr ( void ) const {
	if ( getnrows() != getncols() )
		throw MatrixError();

	Matrix *A = new Matrix ( getnrows(), getncols()*2 );
	Matrix *inv = new Matrix ( getnrows(), getncols() );

	unsigned int i,j,k;

	for ( i=0; i<getnrows(); i++ ) {
		for ( j=0; j<getncols(); j++ ) {
			(*A) (i,j) = (*this)(i,j);
			(*A) (i,j+getncols()) = i==j;
		}
	}

	Matrix * QR = A->qr_dec ();

	//QR->print();

	for ( k=getncols()-1; k<UINT_MAX; k-- ) {
		for ( i=getnrows()-1; i<UINT_MAX; i-- ) {
			for ( j=getncols()-1; j>i; j-- ) {
				(*QR)(i,k+getncols()) -= (*QR)(i,j)*(*QR)(j,k+getncols());
			}
			(*QR)(i,k+getncols()) /= (*QR)(i,i);
			(*inv)(i,k) = (*QR)(i,k+getncols());
		}
		//std::cout << "QR(" << k << ") = ";
		//QR->print();
	}
	delete A;
	delete QR;

	return inv;
}

std::vector<double> Matrix::solve ( const std::vector<double>& b ) {
	Matrix *LU = lu_dec();

	std::vector<double> x ( nrows );
	std::vector<double> y ( nrows );

	y = forward ( LU, b);
	x = backward ( LU, y );

	delete LU;

	return x;
}

std::vector<double> Matrix::forward ( const Matrix *LU, const std::vector<double>& b ) {
	unsigned int i,k;
	double s;
	std::vector<double> y (nrows);

	for (i=0; i<nrows; i++) {
		s = b[i];
		for (k=0; k<i; k++) {
			s -= (*LU)(i,k) * y[k];
		}
		y[i] = s;
	}
	return y;
}

std::vector<double> Matrix::backward ( const Matrix *LU, const std::vector<double>& y ) {
	int i;
	unsigned int k;
	double s;
	std::vector<double> x (nrows);

	for (i=nrows-1; i>=0; i--) {
		s = y[i];
		for (k=i+1; k<nrows; k++) {
			s -= (*LU)(i,k) * x[k];
		}
		x[i] = s/(*LU)(i,i);
	}
	return x;
}

Matrix* Matrix::inverse ( void ) {
	if (nrows!=ncols)
		throw MatrixError();

	Matrix *LU = lu_dec ();
	Matrix *Inv = new Matrix ( nrows, nrows );
	std::vector<double> x ( nrows, 0);
	std::vector<double> y ( nrows, 0);
	unsigned int i,j;

	for ( i=0; i<ncols; i++ ) {
		for (j=0; j<nrows; j++)
			x[j] = 0;
		x[i] = 1;

		y = forward ( LU, x );
		x = backward ( LU, y );

		for (j=0; j<nrows; j++)
			(*Inv)(j,i) = x[j];
	}

	delete LU;

	return Inv;
}

std::vector<double> Matrix::operator* ( std::vector<double>& x ) {
	if (x.size() != ncols)
		throw MatrixError();

	std::vector<double> out ( nrows, 0 );
	unsigned int i,j;

	for (i=0; i<nrows; i++) {
		for (j=0; j<ncols; j++) {
			out[i] += (*this)(i,j) * x[j];
		}
	}

	return out;
}

void Matrix::scale ( double a ) {
	unsigned int i;

	for (i=0; i<nrows*ncols; i++)
		data[i] *= a;
}

bool Matrix::symmetric ( void ) {
	unsigned int i,j;
	for ( i=0; i<nrows; i++ ) {
		for ( j=i; j<ncols; j++ ) {
			if ((*this)(i,j)!=(*this)(j,i))
				return false;
		}
	}
	return true;
}

std::vector<double> leastsq ( const Matrix *A, const std::vector<double>& b ) {
	Matrix *M = new Matrix ( A->getnrows(), A->getncols()+1 );
	unsigned int i, j;
	unsigned int nr (A->getnrows()),nc(A->getncols());

	for ( i=0; i<nr; i++ ) {
		for ( j=0; j<nc; j++ ) {
			(*M)(i,j) = (*A)(i,j);
		}
		(*M)(i,nc) = b[i];
	}

	std::vector<double> out ( leastsq ( M ) );
	delete M;

	return out;
}

std::vector<double> leastsq ( const Matrix *M ) {
	int i,j;
	unsigned int nr(M->getnrows()),nc(M->getncols()-1);
	Matrix *R = M->qr_dec();
	std::vector<double> x (nc);
	double s;

	for ( i=nc-1; i>=0; i-- ) {
		s = (*R)(i,nc);
		for ( j=i+1; j<nc; j++ ) {
			s -= (*R)(i,j)*x[j];
		}
		x[i] = s/(*R)(i,i);
	}

	delete R;

	return x;
}
