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

## MG_DISTS_HYPEREXP
##
## SUMMARY
##  Hyper-Exponential distribution.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

mg_dists_dhyperexp <- function( x, probs, rates, log = FALSE )
{
	if ( length( probs ) != length( rates ) )
	{
		stop( "Probability vector length differ from rates vector length." );
	}

	d <- 0;

	for ( i in 1:length(probs) )
	{
		d <- d + probs[i] * dexp( x, rates[i] );
	}

	if ( log )
	{
		return( log( d ) );
	}

	return( d );
}

mg_dists_phyperexp <- function( q, probs, rates, lower.tail = TRUE, log.p = FALSE )
{
	if ( length( probs ) != length( rates ) )
	{
		stop( "Probability vector length differ from rates vector length." );
	}

	p <- 0;

	for ( i in 1:length(probs) )
	{
		p <- p + probs[i] * pexp( q, rates[i], lower.tail = FALSE );
	}

	if ( lower.tail )
	{
		p <- 1 - p;
	}

	if ( log.p )
	{
		return( log( p ) );
	}

	return( p );
}

mg_dists_qhyperexp <- function ( p, probs, rates, lower.tail = TRUE, log.p = FALSE )
{
	if ( length( probs ) != length( rates ) )
	{
		stop( "Probability vector length differ from rates vector length." );
	}

	# Build the CDF of phase probabilities
	probsCdf <- cumsum( probs );
	# Generates 'length(p)' uniform random numbers (in [0,1])
	u <- runif( length(p) );
	# Finds what probabilities we must choose
	idx <- findInterval( u, probsCdf );
	idx[ idx == 0 ] <- 1;

	# Generates exponential quantiles (each with appropriate rate)
	# and returns them
	return( qexp( p, rates[ idx ], lower.tail = lower.tail, log.p = log.p ) );
}

mg_dists_rhyperexp <- function ( n, probs, rates )
{
	if ( length( probs ) != length( rates ) )
	{
		stop( "Probability vector length differ from rates vector length." );
	}

	# Build the CDF of phase probabilities
	probsCdf <- cumsum( probs );
	# Generates n uniform random numbers (in [0,1])
	u <- runif(n);
	# Finds what probabilities we must choose
	idx <- findInterval( u, probsCdf );

	# Generates n exponential random numbers and returns them
	return( rexp( n, rates[ idx ] ) );
}
