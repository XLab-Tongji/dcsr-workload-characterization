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

## MG_DISTS_STABLE
##
## SUMMARY
##  Stable distributions.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)

## SKEWED STABLE DISTRIBUTIONS ##

source("lib/mg/mg_dists_stable_fBasics_2A_StableDistribution.R");
source("lib/mg/mg_dists_stable_fBasics_2D_DistributionFits.R");

mg_dists_dstable <- function( x, alpha, beta, gamma = 1, delta = 0, pm = c(0, 1, 2), log = FALSE )
{
#	havefBasics <- require( "fBasics" );
#
#	if ( havefBasics )
#	{
#		d <- fBasics::dstable( x, alpha, beta, gamma, delta, pm );
#		if ( log )
#		{
#			d <- log( d );
#		}
#		return( d );
#	}
#
#	stop( "Not yet implemented!" );

	d <- dstable( x, alpha, beta, gamma, delta, pm );
	if ( log )
	{
		d <- log( d );
	}
	return( d );
}

mg_dists_pstable <- function( q, alpha, beta, gamma = 1, delta = 0, pm = c(0, 1, 2), lower.tail = TRUE, log.p = FALSE )
{
#	havefBasics <- require( "fBasics" );
#
#	if ( havefBasics )
#	{
#		p <- fBasics::pstable( q, alpha, beta, gamma, delta, pm );
#		if ( !lower.tail )
#		{
#			p <- 1 - p;
#		}
#		if ( log.p )
#		{
#			p <- log( p );
#		}
#		return( p );
#	}
#
#	stop( "Not yet implemented!" );

	p <- pstable( q, alpha, beta, gamma, delta, pm );
	if ( !lower.tail )
	{
		p <- 1 - p;
	}
	if ( log.p )
	{
		p <- log( p );
	}
	return( p );
}

mg_dists_qstable <- function( p, alpha, beta, gamma = 1, delta = 0, pm = c(0, 1, 2), lower.tail = TRUE, log.p = FALSE )
{
#	havefBasics <- require( "fBasics" );
#
#	if ( havefBasics )
#	{
#		if ( log.p )
#		{
#			p <- exp( p );
#		}
#		if ( !lower.tail )
#		{
#			p <- 1 - p;
#		}
#		return( fBasics::qstable( p, alpha, beta, gamma, delta, pm ) );
#	}
#
#	stop( "Not yet implemented!" );

	if ( log.p )
	{
		p <- exp( p );
	}
	if ( !lower.tail )
	{
		p <- 1 - p;
	}
	return( qstable( p, alpha, beta, gamma, delta, pm ) );
}

mg_dists_rstable <- function( n, alpha, beta, gamma = 1, delta = 0, pm = c(0, 1, 2) )
{
#	havefBasics <- require( "fBasics" );
#
#	if ( havefBasics )
#	{
#		return( fBasics::rstable( n, alpha, beta, gamma, delta, pm ) );
#	}
#
#	stop( "Not yet implemented!" );

	return( rstable( n, alpha, beta, gamma, delta, pm ) );
}

## SYMMETRIC STABLE DISTRIBUTIONS ##

mg_dists_dsymstb <- function( x, alpha = 1.8, log = FALSE )
{
#	havefBasics <- require( "fBasics" );
#
#	if ( havefBasics )
#	{
#		d <- fBasics::dsymstb( x, alpha );
#		if ( log )
#		{
#			d <- log( d );
#		}
#		return( d );
#	}
#
#	stop( "Not yet implemented!" );

	d <- dsymstb( x, alpha );
	if ( log )
	{
		d <- log( d );
	}
	return( d );
}

mg_dists_psymstb <- function( q, alpha = 1.8, lower.tail = TRUE, log.p = FALSE )
{
#	havefBasics <- require( "fBasics" );
#
#	if ( havefBasics )
#	{
#		p <- fBasics::psymstb( q, alpha );
#		if ( !lower.tail )
#		{
#			p <- 1 - p;
#		}
#		if ( log.p )
#		{
#			p <- log( p );
#		}
#		return( p );
#	}
#
#	stop( "Not yet implemented!" );

	p <- psymstb( q, alpha );
	if ( !lower.tail )
	{
		p <- 1 - p;
	}
	if ( log.p )
	{
		p <- log( p );
	}
	return( p );
}

mg_dists_qsymstb <- function( p, alpha = 1.8, lower.tail = TRUE, log.p = FALSE )
{
#	havefBasics <- require( "fBasics" );
#
#	if ( havefBasics )
#	{
#		if ( log.p )
#		{
#			p <- exp( p );
#		}
#		if ( !lower.tail )
#		{
#			p <- 1 - p;
#		}
#		return( fBasics::qsymstb( p, alpha ) );
#	}
#
#	stop( "Not yet implemented!" );

	if ( log.p )
	{
		p <- exp( p );
	}
	if ( !lower.tail )
	{
		p <- 1 - p;
	}
	return( qsymstb( p, alpha ) );
}

mg_dists_rsymstb <- function( n, alpha = 1.8 )
{
#	havefBasics <- require( "fBasics" );
#
#	if ( havefBasics )
#	{
#		return( fBasics::rsymstb( n, alpha ) );
#	}
#
#	stop( "Not yet implemented!" );

	return( rsymstb( n, alpha ) );
}

mg_dists_stable.fit.mle <- function( x, doplot = FALSE )
{
	fit <- stableFit( x, type = c( "mle" ), doplot = doplot );

	return(
		list(
			alpha = fit@fit$estimate["alpha"],
			beta = fit@fit$estimate["beta"],
			gamma = fit@fit$estimate["gamma"],
			delta = fit@fit$estimate["delta"]
		)
	);
}

mg_dists_stable.fit.quantile <- function( x, doplot = FALSE )
{
	fit <- stableFit( x, type = c( "q" ), doplot = doplot );

	return(
		list(
			alpha = fit@fit$estimate["alpha"],
			beta = fit@fit$estimate["beta"],
			gamma = fit@fit$estimate["gamma"],
			delta = fit@fit$estimate["delta"]
		)
	);
}
