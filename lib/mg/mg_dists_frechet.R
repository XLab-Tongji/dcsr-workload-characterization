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

## MG_DISTS_FRECHET
##
## SUMMARY
##  Frechet (or Type II Extreme Value) distribution.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

mg_dists_dfrechet <- function( x, location = 0, scale = 1, shape = 1, log = FALSE )
{
	mg_dists_frechet_haveevd <- require( "evd" );

	if ( mg_dists_frechet_haveevd )
	{
		return( evd::dfrechet( x, location, scale, shape, log ) );
	}

	# preconditions
	if ( scale <= 0 )
	{
		stop( "[MG::Dists::Frechet::DFrechet] Scale parameter is out of range." );
	}
	if ( shape <= 0 )
	{
		stop( "[MG::Dists::Frechet::DFrechet] Shape parameter is out of range." );
	}
	if ( all( x <= (location - scale * shape) ) )
	{
		stop( "[MG::Dists::Frechet::DFrechet] Data are out of range." );
	}

	ix = which( x > 0, arr.ind = TRUE );

	z <- (x[ix] - location) / scale;
	d[ix] <- log(shape / scale) - (1 + shape) * log( z ) - z^(-shape);

	if ( !log )
	{
		d[ix] <- exp( d[ix] );
	}

	d[ x<= 0 ] <- -Inf;

	return( d );
}

mg_dists_pfrechet <- function( q, location = 0, scale = 1, shape = 1, lower.tail = TRUE, log.p = FALSE )
{
	mg_dists_frechet_haveevd <- require( "evd" );

	if ( mg_dists_frechet_haveevd )
	{
		p <- evd::pfrechet( q, location, scale, shape, lower.tail );
		if ( log.p )
		{
			p <- log( p );
		}
		return( p );
	}

	stop( "Not yet implemented!" ); #FIXME
}

mg_dists_qfrechet <- function( p, location = 0, scale = 1, shape = 1, lower.tail = TRUE, log.p = FALSE )
{
	mg_dists_frechet_haveevd <- require( "evd" );

	if ( mg_dists_frechet_haveevd )
	{
		if ( log.p )
		{
			p <- exp( p );
		}
		return( evd::qfrechet( p, location, scale, shape, lower.tail ) );
	}

	stop( "Not yet implemented!" ); #FIXME
}

mg_dists_rfrechet <- function( n, location = 0, scale = 1, shape = 1 )
{
	mg_dists_frechet_haveevd <- require( "evd" );

	if ( mg_dists_frechet_haveevd )
	{
		return( evd::rfrechet( n, location, scale, shape ) );
	}

	if ( any( scale <= 0 ) )
	{
		stop( "[MG::Dists::Frechet::RFrechet] Scale parameter is out-of-range (<=0)." );
	}
	if ( any( shape <= 0 ) )
	{
		stop( "[MG::Dists::Frechet::RFrechet] Shape parameter is out-of-range (<=0)." );
	}

	return( location + scale * (- log( runif(n) ) )^(- shape ) );
}

## Fits a Frechet distribution with Maximum Likelihood method
mg_dists_frechet.fit.mle <- function( x )
{
	n <- length(x);

	# Initial values
	# FIXME: how to choose good initial values?
	#sigma0 <- ( mu0 - mean(x) ) / (1 - gamma(1 - 1/alpha0)); # Derived from E[X] expression (we can use it 'cause alpha0 > 1)
	#sigma0 <- 1;
	sigma0 <- sqrt(6) * sd(x) / pi;
	#mu0 <- 0;
	mu0 <- mean(x) - 0.57722 * sigma0;

	#alpha0 <- .1; # alpha > 2 => variance is finite
	alpha0 <- 2 + .Machine$double.eps^.25; # alpha > 2 => variance is finite

	parms0 <- c( sigma0, mu0, alpha0 );

	# The log-likelihood function
	negloglik <- function( parmshat, x )
	{
		# parmshat == c( sigma, mu, alpha )

		if ( parmshat[1] <= 0 || parmshat[3] <= 0 )
		{ 
			res <- 1e+06;
		} else {
			y <- (x - parmshat[2]) / parmshat[1];

			term1 <- length(x) * log( parmshat[3]/parmshat[1] );
			term2 <- (1 + parmshat[3]) * sum( log( y ) );
			term3 <- sum( y^(- parmshat[3] ) );
			res <- -term1 + term2 + term3;
		}
		return( res );
	};

	# Solve log-Likelihood(parmshat) = 0
	fit <- optim( par = parms0, negloglik, x = x );
	if ( fit$convergence )
	{
		warning( "[MG::Dists::Frechet::Fit::MLE] Optimization may not have succeeded." );
	}
	return( list( location = fit$par[2], scale = fit$par[1], shape = fit$par[3]  ) );
}
