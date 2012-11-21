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
##  A collection of vector utilities.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

## Returns a vector where all elements are 1.
mg_vector_ones <- function( n )
{
	return( array( 1, n ) );
}

## Returns a vector where all elements are 0.
mg_vector_zeros <- function( n )
{
	return( array( 0, n ) );
}

mg_vector_t <- function( v )
{
	v <- as.vector(v);
	if ( is.null( dim(v) ) )
	{
		# v is a row vector
		return( as.matrix(v) );
	}

	# v is a column vector
	return( as.vector( v ) );
}

## Returns the "dot" (or "scalar") product "X . Y", where "X" and "Y" are
## vectors of any length.
mg_vector_dot <- function( x, y )
{
	#return( sum( x * y ) );
	#return( t(x) %*% y );
	#return( crossprod( x, y ) );
	return( as.vector(x) %*% as.vector(y) ); # same as above
}

## Returns the "cross" (or "vector") product "X x Y", where "X" and "Y" are
## vectors of size 2 or 3.
mg_vector_cross <- function( x, y )
{
	x <- as.vector(x);
	y <- as.vector(y);
	nx <- length(x);
	ny <- length(y);

	# preconditions
	if ( (nx != 2 && nx != 3) && (ny != 2 && ny != 3) )
	{
		stop( "[MG::Vector::Cross] Dimension of input vectors must be 2 or 3." );
	}

	res <- NULL;
	if ( nx == 2 )
	{
		if ( ny == 2 )
		{
			res <- c( x[1]*y[2] - x[2]*y[1] );
		} else {
			res <- c( x[2]*y[3], -x[1]*y[3], x[1]*y[2] - x[2]*y[1] );
		}
	} else {
		if ( ny == 2 )
		{
			res <- c( -x[3]*y[2], x[3]*y[1], x[1]*y[2] - x[2]*y[1] );
		} else {
			res <- c( x[2]*y[3] - x[3]*y[2], x[3]*y[1] - x[1]*y[3], x[1]*y[2] - x[2]*y[1] );
		}
	}

	return( res );
}

## Retuns "a * X", where "a" is a scalar and "X" is a vector.
mg_vector_scale <- function( x, a )
{
	return( as.numeric(a) * as.vector(x) );
}
