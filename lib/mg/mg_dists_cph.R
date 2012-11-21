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

## MG_DISTS_CPH
##
## Continuos Phase-Type PH(alpha, T) distribution.
## 
## * alpha: initial probability distribution (row vector of length n)
## * T: generator matrix (nxn matrix)
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

source( "lib/mg/mg_matrix.R" );
source( "lib/mg/mg_rvg.R" );
source( "lib/mg/mg_vector.R" );

##
## Returns the value of density function:
##  f(x) = alpha*exp{T*x}*t
## where:
## * "t=-T*ones(n,1)"
## * "ones(n,1)" is a column vector of length "n" of ones
## * "n" is the length of "alpha" vector.
##
mg_dists_dcph <- function( x, alpha, T, log = FALSE )
{
	if ( length( alpha ) != ncol( T ) )
	{
		stop( "The size of 'alpha' vector must equal the number of columns of 'T' matrix." );
	}
	if ( ncol(T) != nrow(T) )
	{
		stop( "The 'T' matrix must be squared." );
	}

	#ones <- matrix( 1, nrow = length( alpha ), ncol = 1 ); # a column vector of 1s
	ones <- mg_vector_ones( length( alpha ) ); # a column vector of 1s
	#t <- - T %*% ones; # the exit vector
	t <- - mg_matrix_prodv( T, ones ); # the exit vector

#	# coercion from "matrix" to "Matrix"
#	if ( !is( T, "Matrix" ) )
#	{
#		T <- as( T, "Matrix" );
#	}

	# trasform initial vector to a row vector
	#alpha <-  matrix(alpha, nrow = 1, ncol = length(alpha) );

	d <- c();
	for ( i in 1:length(x) )
	{
		#d <- append( d, as.numeric( alpha %*% Matrix::expm( T * x[i] ) %*% t ) );
		d <- append( d, mg_vector_dot( mg_matrix_vprod( alpha, mg_matrix_exp( T * x[i] ) ), t ) );
	}

	if ( log )
	{
		d <- log( d );
	}

	return( d );
}

##
## Returns the value of cumulative distribution function:
##  F(x) = P(X <= x) = 1-alpha*exp{T*x}*ones(n,1)
## where:
## * "ones(n,1)" is a column vector of length "n" of ones
## * "n" is the length of "alpha" vector.
##
mg_dists_pcph <- function( q, alpha, T, lower.tail = TRUE, log.p = FALSE )
{
	if ( length( alpha ) != ncol( T ) )
	{
		stop( "The size of 'alpha' vector must equal the number of columns of 'T' matrix." );
	}

	#ones <- matrix( 1, nrow = length( alpha ), ncol = 1 ); # a column vector of 1s
	ones <- mg_vector_ones( length(alpha) );
	#t <- - T %*% ones; # the exit vector
	t <- - mg_matrix_prodv( T, ones ); # the exit vector

#	# coercion from "matrix" to "Matrix"
#	if ( !is( T, "Matrix" ) )
#	{
#		T <- Matrix( T );
#	}

	# trasform initial vector to a row vector
	#alpha <-  matrix(alpha, nrow = 1, ncol = length(alpha) );

	p <- c();
	for ( i in 1:length(q) )
	{
		##p <- append( p, as.numeric( alpha %*% Matrix::expm( T * q[i] ) %*% ones ) );
		#p <- append( p, mg_vector_dot( mg_matrix_vprod( alpha, Matrix::expm( T * q[i] ) ), ones ) );
		p <- append( p, mg_vector_dot( mg_matrix_vprod( alpha, mg_matrix_exp( T * q[i] ) ), ones ) );
	}

	if ( lower.tail )
	{
		p <- 1 - p;
	}

	if ( log.p )
	{
		p <- log( p );
	}

	return( p );
}

##
## Returns the quantile value:
##  q = F^{-1}(p)
## where:
## * "F^{-1}(p)" is the inverse cumulative distribution function
##
mg_dists_qcph <- function( p, alpha, T, lower.tail = TRUE, log.p = FALSE )
{
	if ( log.p )
	{
		p <- exp(p);
	}
	if ( lower.tail )
	{
		p <- 1-p;
	}

	N_IT <- 100; # max iteration if tolerance is not reached
	#N_IT <- 1; # max iteration if tolerance is not reached
	n <- length(p);
	x <- mg_dists_cph.moment( 1, alpha, T );
	x <- rep( x, n ); # transform x to a vector
	z <- mg_vector_zeros( n );
	#tol <- 1.0e-10;
	tol <- .Machine$double.eps^.5; # a reasonable accurary for equality
	#tol <- .Machine$double.eps;
	for ( i in 1:N_IT )
	{
		derv <- mg_dists_dcph( x, alpha, T );
		dif <- mg_dists_pcph(x, alpha, T) - p;

		#if( all( abs( dif ) <= tol ) )
		if( all( abs( dif[ !is.na(dif) ] ) <= tol ) )
		{
			return( x );
		}
		x <- pmax( x - dif / derv, z );
	}

	warning( "[MG::DISTS::CPH::QCPH] Accuracy not reached, returns zeros" );

	return( z );
}

##
## Returns "n" random numbers generated from the given "PH(alpha,T)"
## distribution.
##
mg_dists_rcph <- function( n, alpha, T, ... )
{
	return( mg_dists_rcph._neutspagano( n, alpha, T, ... )$r );
}

##
## Returns "n" random numbers generated from the given "PH(alpha,T)"
## distribution.
##
## Same as "mg_dists_rcph" but also returns the random variate generator (rvg)
## internal state.
##
mg_dists_rcph.ex <- function( n, alpha, T, ... )
{
	return( mg_dists_rcph._neutspagano( n, alpha, T, ... ) );
}

##
## Returns the non-central kth-moment:
##  E[X^k] = (-1)^k * k! * alpha * T^(-k) * I
##
## TODO: We found two possible methods. Test what of two is better (i.e.
##       faster, more stable, use less memory).
##
mg_dists_cph.moment <- function( k, alpha, T )
{
	n <- length( alpha );

#BEGIN: First method
	return(
		as.numeric(
			(-1)^k
			* factorial(k)
			* mg_vector_dot(
				mg_matrix_vprod(
					alpha,
					mg_matrix_pow(
						mg_matrix_inv( T ),
						k
					)
				),
				mg_vector_ones(n)
			)
		)
	);
#END: First method

#BEGIN: Second method
#	res <- 0;
#	x <- mg_vector_zeros( n );
#	T1 <- mg_matrix_t( mg_matrix_pow( T, k ) );
#	TT <- mg_matrix_prod( mg_matrix_t( T1 ), T1 );
#	b <- mg_matrix_tprodv( T1, alpha );
#	x <- solve( TT, b );
#	res <- mg_vector_dot( x, mg_vector_ones( n ) );
#	if ( (k %% 2) != 0 )
#	{
#		res <- -res;
#	}
#	res <- res * factorial(k);
#	return( as.numeric( res ) );
#END: Second method
}

##
## Returns a random numbers according to the given "PH(alpha,T)"
## distribution, generated using the Neuts-Pagano method.
##
## SEE ALSO:
##  * M. Neuts, M. E. Pagano. "Generating Random Variates from a Distribution of Phase Type". (1981)
##
mg_dists_rcph._neutspagano <- function( n, alpha, T, rvg.state )
{
	if ( missing( rvg.state ) )
	{
		rvg.state <- mg_dists_rcph._neutspagano.init( alpha, T );
	}

	m <- length( alpha );
	alphax <- c( alpha, 1 - sum( alpha ) ); # alpha + absorbing state

	h <- 0;
	xph <- c();
	while ( h <= n )
	{
		# Simulates the underlying Markov chain

		xph[h] <- 0;
		#k <- array( 0, m );
		k <- rep.int( 0, m );

		# From (alpha, absorbing_state) choose current state i
		i <- mg_rvg_alias.rand( 1, alphax, rvg.state$alpha.aliascut );

		absorbed <- FALSE;
		while ( !absorbed )
		{
			if ( i == (m+1) )
			{
				absorbed <- TRUE;
				break;
			}
			k[i] <- k[i] + 1;

			# From hatQ[i] choose current state j
			j <- mg_rvg_alias.rand( 1, rvg.state$hatQ[i,], rvg.state$hatQ.aliascut[[i]] );

			i <- j;
		}

		for ( i in 1:m )
		{
			if ( k[i] > 0 )
			{
				# state i was visited at least one time

				xph[h] <- xph[h] + mg_dists_cph.rerlang( 1, k[i], rvg.state$xexp[i] );
				k[i] <- 0;
			}
		}

		h <- h + 1;
	}

	#return( xph );
	return( list( r = xph, rvg.state = rvg.state ) );
}

##
## Initialize the alias-cutoff tables for the Neuts-Pagano method.
##
mg_dists_rcph._neutspagano.init <- function( alpha, T )
{
	m <- length( alpha );
	alphax <- c( alpha, 1 - sum( alpha ) ); # alpha + absorbing state

	# Computes alpha aliat-cutoff info
	alpha.ac <- mg_rvg_alias.init( alphax );

	# Computes hatQ and haQ alias-cutoff info
	xexp <- c(); # stores -Q[i,i]^-1; used for generating Erlang variates
	hatQ <- matrix( 0, m, m+1 );
	hatQ.ac <- list();
	for ( i in 1:m )
	{
		# Computes hatQ[i,.]
		xexp[i] <- -T[i,i];
		for ( j in 1:m )
		{
			hatQ[i,j] <- T[i,j] * xexp[i];
		}
		hatQ[i,i] <- 0;
		hatQ[i,m+1] <- 1 - min( sum( hatQ[i,] ) );

		ac <- mg_rvg_alias.init( hatQ[i,] );
		hatQ.ac[[i]] <- ac;
	}

	return( list( xexp = xexp, hatQ = hatQ, alpha.aliascut = alpha.ac, hatQ.aliascut = hatQ.ac ) );
}

##
## Generates "n" random numbers according to the given
## "Erlang(k,lambda)" distribution.
##
mg_dists_cph.rerlang <- function( n, k, lambda )
{
	if ( k <= 0 )
	{
		warning( "[MG::DISTS::CPH::RCPH::RERLANG] FIXME k <= 0" ); #FIXME
		return( 0 );
	}

	r <- c();
	for ( j in 1:n  )
	{
		u <- 1;
		for ( i in 1:k )
		{
			u <- u * runif(1);
		}

		r[j] <- - (1/lambda) * log(u);
	}

	return( r );
}
