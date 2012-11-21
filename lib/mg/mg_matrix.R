##
## Copyright (C) 2008       Marco Guazzone
##                          [Distributed Computing System (DCS) Group,
##                           Computer Science Institute,
##                           Department of Science and Technological Innovation,
##                           University of Piemonte Orientale,
##                           Alessandria (Italy)]
##
## This file is part of dcsr-workload-characterization
##
## dcsr-workload-characterization is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## dcsr-workload-characterization is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with dcsr-workload-characterization  If not, see <http://www.gnu.org/licenses/>.
##

## MG_MATRIX
##
## SUMMARY
##  A collection of matrix utilities.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

source( "lib/mg/mg_math_funcs.R" );
source( "lib/mg/mg_vector.R" );

## Returns a "nrows x ncols" matrix of 1.
mg_matrix_ones <- function( nrows, ncols )
{
	return( matrix( 1, nrow = nrows, ncol = ncols ) );
}

## Returns a "nrows x ncols" matrix of 0.
mg_matrix_zeros <- function( nrows, ncols )
{
	return( matrix( 0, nrow = nrows, ncol = ncols ) );
}

## Returns a vector containing the number of row and column, respectively, if
## parameter "d" is not specified. Otherwise, returns the size of the specified
## dimension.
mg_matrix_size <- function( A, d = NULL )
{
	if ( is.null( d ) )
	{
		return( dim( as.matrix(A) ) );
	}

	return( dim( as.matrix(A) )[ d ] );
}

## Returns the "A'" (i.e. the transpose of "A"), where "A" is a matrix.
mg_matrix_t <- function( A )
{
	return( t( as.matrix(A) ) );
}

## Returns "A^(-1)", i.e. the inverse of "A", where "A" is a matrix.
mg_matrix_inv <- function( A )
{
	return( solve( as.matrix(A) ) );
}

## Returns "a x A", where "a" is a scalar and "A" is a matrix.
mg_matrix_scale <- function( a, A )
{
	return( as.numeric(a) * as.matrix(A) );
}

## Returns "A^k", where "A" is a matrix.
mg_matrix_pow <- function( A, k )
{
	return( A^k );
}

## Returns "A x B", where "A" and "B" are matrices of proper size.
mg_matrix_prod <- function( A, B )
{
	return( A %*% B );
}

## Returns "A' x B", where "A" and "B" are matrix of proper size.
mg_matrix_tprod <- function( A, B )
{
	return( mg_matrix_prod( t(A), B ) );
}

## Returns "A x v", where "A" is matrix and "v" is a vector of proper size.
mg_matrix_prodv <- function( A, v )
{
	return( A %*% v );
}

## Returns "A' x v", where "A" is matrix and "v" is a vector of proper size.
mg_matrix_tprodv <- function( A, v )
{
	return( mg_matrix_prodv( t(A), v ) );
}

## Returns "v x A", where "A" is matrix and "v" is a vector of proper size.
mg_matrix_vprod <- function( v, A )
{
	return( as.vector(v) %*% as.matrix(A) );
}

## Returns "v x A'", where "A" is matrix and "v" is a vector of proper size.
mg_matrix_vprodt <- function( v, A )
{
	return( mg_matrix_prodv( v, as.matrix(t(A)) ) );
}

## Returns "v' x A", where "A" is matrix and "v" is a vector of proper size.
mg_matrix_vtprod <- function( v, A )
{
	return( mg_vector_t(v) %*% as.matrix(A) );
}

##
## Returns the matrix exponential:
##  exp( A ) = e^A
## i.e. the constant "e" to the matrix power "A".
##
mg_matrix_exp <- function( A )
{
	#return( mg_matrix_exp.taylor( A ) );
	#eturn( mg_matrix_exp.eigenv( A ) );
	#return( mg_matrix_exp.pade.Matrix( A ) );
	return( mg_matrix_exp.pade.C( A ) );
	#return( mg_matrix_exp.pade( A ) );
}

##
## Returns the matrix exponential via Taulor series.
##
## As a practical numerical method, this is often slow and inaccurate.
mg_matrix_exp.taylor <- function( A )
{
	haveMatrix <- require( Matrix );

	if ( !haveMatrix )
	{
		stop( "Package 'Matrix' required." );
	}

	# coercion from "matrix" to "Matrix"
	if ( !is( A, "Matrix" ) )
	{
		A <- as( A, "Matrix" );
	}

	return( Matrix::expm( A ) );
}

##
## Returns the matrix exponential via Pade' approximation with scaling and
## squaring. 
##
## SEE ALSO
## * Golub, G. H. and C. F. Van Loan, Matrix Computation, p. 384, Johns Hopkins University Press, 1983.
## * Moler, C. B. and C. F. Van Loan, "Nineteen Dubious Ways to Compute the Exponential of a Matrix," SIAM Review 20, 1978, pp. 801-836.
## * Higham, N. J., "The Scaling and Squaring Method for the Matrix Exponential Revisited," SIAM J. Matrix Anal. Appl., 26(4) (2005), pp. 1179-1193.
##
mg_matrix_exp.pade <- function( A )
{
	# Scale A by power of 2 so that its norm is < 1/2 .
	lfe <- mg_math_log2( mg_matrix_norm(A, Inf) );
	s <- max(0,lfe$e+1);
	A <- A/2^s;

	# Pade approximation for exp(A)
	X <- A;
	c <- 1/2;
	E <- diag(1, nrow(A), ncol(A)) + c*A;
	D <- diag(1, nrow(A), ncol(A)) - c*A;
	q <- 6;
	p <- TRUE;
	for ( k in 2:q )
	{
		c <- c * (q-k+1) / (k*(2*q-k+1));
		X <- A %*% X;
		cX <- c * X;
		E <- E + cX;
		if ( p )	
		{
			D <- D + cX;
		} else {
			D <- D - cX;
		}
		p <- !p;
	}
	E <- solve(D, E);

	# Undo scaling by repeated squaring
	for ( k in 1:s )
	{
		E <- E %*% E;
	}

	return( E );
}

## C implementation (it seems slower than R's version)
mg_matrix_exp.pade.C <- function( A )
{
	dyn.load( "lib/mg/src/mg_matrix.so" );

	return( .Call( "mg_matrix_exp", A ) );
}

mg_matrix_exp.pade.Matrix <- function( A )
{
	haveMatrix <- require( Matrix );

	if ( !haveMatrix )
	{
		stop( "[MG::MATRIX::EXP::PADE::MATRIX] Package 'Matrix' required." );
	}

	# coercion from "matrix" to "Matrix"
	if ( !is( A, "Matrix" ) )
	{
		A <- as( A, "dgeMatrix" );
	}

	return( as.matrix( Matrix::expm( A ) ) );
}

##
## Returns the matrix exponential via eigenvalues and eigenvectors.
##
## Method:
##   exp(A) = V*exp(D)*V^(-1)
## where:
##   * A is a square matrix of order n
##   * V is a matrix whose column are eigenvectors of A
##   * D is a diagonal matrix diag{d_1,...,d_n} where d_i is a eigenvalue of A
##     (and v_{.j} is a related eigenvector).
##
## As a practical numerical method, the accuracy is determined by the condition
## of the eigenvector matrix.
##
mg_matrix_exp.eigenv <- function( A )
{
	e <- eigen( A );
	return( e$vectors %*% diag( exp( e$values ) ) %*% solve( e$vectors ) );
}

## MG_MATRIX_LOG.EIGENV
##
## SUMMARY
##  Matrix logarithm (eigenvalue version).
##
mg_matrix_log.eigenv <- function( A )
{
	e <- eigen( A );
	return( e$vectors %*% diag( log( diag( e$values ) ) ) %*% solve( e$vectors ) );
}

##
## Computes the p-norm of the matrix, or the vector, @c A.
##
## For matrices...
##  * @c p = 2:     is the 2-norm of A, largest singular value of A, max(svd(A)).
##  * @c p = 1:     is the 1-norm of A, the largest column sum, max(sum(abs(A))).
##  * @c p = Inf:   is the infinity norm of A, the largest row sum, max(sum(abs(t(A)))).
##  * @c p = 'fro': is the Frobenius norm, sqrt(sum(diag(t(A)%*%A))).
##  * @c p = N:      is available for matrix A only if N is 1, 2, Inf or 'fro'.
##
## For vectors...
##  * @c p = N:    sum(abs(A)^N)^(1/N).
##  * @c p = Inf:  max(abs(A)).
##  * @c p = -Inf: min(abs(A)).
##
mg_matrix_norm <- function( A, p = 2 )
{
	if ( is.vector( A ) )
	{
		if ( p == "fro" )
		{
			return( sqrt( sum( abs( A )^2 ) ) );
		} else if ( p == Inf ) {
			return( max( abs( A ) ) );
		} else if ( p == -Inf ) {
			return( min( abs( A ) ) );
		} else if ( p != 0 ) {
			return( sum( abs( A )^p )^(1/p) );
		} else {
			stop( "[MG::Matrix::Norm] Unrecognized p-norm, with p=0." );
		}
	} else if ( is.matrix( A ) ) {
		if ( p == "fro" )
		{
			#return( sqrt( max( sum( diag( t(A) %*% A ) ) ) ) );
			return( sqrt( sum( colSums( abs( A )^2 ) ) ) );
		}
		if ( p == Inf )
		{
			return( max( colSums( abs( t(A) ) ) ) );
		}
		if ( p == 1 )
		{
			return( max( colSums( abs( A ) ) ) );
		}
		if ( p == 2 )
		{
			return( max( svd( A )$d ) );
		}
	}

	stop( "[MG::Matrix::Norm] Unrecognized p-norm." );
}
