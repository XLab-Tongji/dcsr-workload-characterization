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

## MG_DISTS_PARETO
##
## SUMMARY
##  Pareto distribution.
##
## DESCRIPTION
##  Pareto distribution:
##   f(x)=\frac{k\,x_\mathrm{m}^k}{x^{k+1}}
##   F(x)=1-\left(\frac{x_\mathrm{m}}{x}\right)^k
##  Parameters:
##  * k: the shape parameter.
##  * x_m: the minimum value (also called location).
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

mg_dists_dpareto <- function( x, shape, xmin = .Machine$double.eps^.25, log = FALSE )
{
	if ( shape <= 0 )
	{
		return( rep( NaN, length(x) ) );
	}

	# out-of-range values will get a NaN
	d <- rep( NaN, length(x) );

	# Find in-range values
	good <- x >= xmin | is.na(x);
	#d[!good] <- NaN;

	#d <- shape * ( xmin^shape ) / ( x^( shape + 1) );
	d[good] <- log(shape) + shape*log(xmin) - (shape + 1)*log(x[good]); # log( f(x) )

	if ( !log )
	{
		d[good] <- exp( d[good] );
	}

	return( d );
}

mg_dists_ppareto <- function( q, shape, xmin = .Machine$double.eps^.25, lower.tail = TRUE, log.p = FALSE )
{
	if ( shape <= 0 )
	{
		return( rep( NaN, length(q) ) );
	}

	# out-of-range values will get a zero
	p <- rep( 1, length(q) );

	# Find in-range values
	good <- q >= xmin | is.na(q);
	#p[!good] <- 0;

	#p <- ( xmin / q )^shape; # 1 - F(q)
	p[good] <- shape * (log(xmin) - log(q[good])); # log( 1 - F(q) )

	if ( !log.p )
	{
		# log( 1-F(q) ) => 1-F(q)
		p[good] <- exp( p[good] );
	}

	if ( lower.tail )
	{
		#p <- 1 - ( xmin / q )^shape;
		p[good] <- 1 - p[good]; # F(x)
	}
#	else
#	{
#		p <- ( xmin / q )^shape;
#	}

	return( p );
}

mg_dists_qpareto <- function ( p, shape, xmin = .Machine$double.eps^.25, lower.tail = TRUE, log.p = FALSE )
{
	if ( shape <= 0 )
	{
		return( rep( NaN, length(p) ) );
	}

	if ( log.p )
	{
		p <- exp( p );
	}

	# out-of-range values will get a NaN
	q <- rep( NaN, length(p) );

	good <- (p >= 0 && p <= 1) | is.na(p);
	#q[!good] <- NaN;

	if ( !lower.tail )
	{
		p[good] <- 1 - p[good];
	}

	q[good] <- xmin / ( ( 1 - p[good] )^( 1 / shape ) );

	return( q );
}

mg_dists_rpareto <- function ( n, shape, xmin = .Machine$double.eps^.25 )
{
	if ( shape <= 0 )
	{
		return( rep( NaN, n ) );
	}

	return( mg_dists_qpareto( runif(n), shape = shape, xmin = xmin ) );
}
