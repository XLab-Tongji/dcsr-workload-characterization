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

## MG_UTILS
##
## SUMMARY
##   Generic utility functions.
##
## DESCRIPTION
##   A set of useful functions.
##
## AUTHORS
##   Marco Guazzone (marco.guazzone@gmail.com)
##

## MG_UTILS_RMALL
##
## Removes all objects
##
## FIXME: it seems this function doesn't work; however cutting&pasting its
##        instructions on R prompt, it works! Why?
mg_utils_rmall <- function()
{
	rm( list = ls() );
}

## MG_UTILS_CLONE
##
## Clones an object n times and returns all occurrences in a vector.
mg_utils_clone <- function( x, n = 1 )
{
	#seq( from = x, to = x, length.out = n );
	return( rep( x, n ) );
}

## MG_UTILS_SEARCH.PATH
##
# search for file in paths
# fn is filename
# paths is a vector of path names, default is constructed from PATH
#
# e.g.: source(search.path("myscript.R"))
mg_utils_search.path <- function( fn, paths )
{
	if ( missing( paths ) )
	{
		paths = strsplit(
			chartr( "\\", "/", Sys.getenv("PATH") ),
			split = switch( .Platform$OS.type, windows = ";", ":")
		)[[1]];
	}

	for ( d in paths )
	{
		if ( file.exists( f <- file.path(d, fn) ) )
		{
			return( f );
		}
	}
	return( NULL );
}

mg_utils_lookup <- function( x, y, leftmost.closed = FALSE, rightmost.closed = FALSE )
{
	res <- c();

	res <- findInterval( y, x );
	if ( leftmost.closed )
	{
		# Replaces 0 with 1
		res <- pmax( res, 1 );
	}
	if ( rightmost.closed )
	{
		# Replaces n+1 with n
		res <- pmin( res, length(y) );
	}

	return( res );
}

mg_utils_lookup.old <- function( x, y )
{
	idx <- c();

	if ( is.vector( x ) )
	{
		ydim <- c(1, 1);
		if ( is.vector( y ) )
		{
			ydim[1] <- length( y );
		} else if ( is.matrix( y ) ) {
			ydim <- dim( y );
		}

		xlen <- length( x );

		if ( x[1] > x[ xlen ] )
		{
			## decreasing table
			p <- order( c( as.vector( y ), as.vector( x ) ) );
			idx[ p ] <- cumsum( as.numeric( p > ydim[1] * ydim[2] ) );
			idx <- xlen - idx[ 1 : ( ydim[1]*ydim[2] ) ];
		} else {
			## increasing table
			p <- order( c( as.vector( x ), as.vector( y ) ) );
			idx[ p ] <- cumsum( as.numeric( p <= xlen ) );
			idx <- idx[ ( xlen + 1 ) : ( xlen + ydim[1]*ydim[2] ) ];
		} 
		#idx <- matrix( idx, nrow = ydim[1], ncol = ydim[2], byrow = TRUE );
	}

	return( idx );
}

## data <- mg_utils_sort.data.frame(data, key = "mg")
## TODO: treat NA (e.g. na.last = TRUE|FALSE|NA
mg_utils_sort.data.frame <- function(x, key, ...) {
	if ( missing(key) )
	{
		rn <- rownames(x);
		if ( all( rn %in% 1:nrow(x) ) )
		{
			rn <- as.numeric(rn)
		}
		x[ order( rn, ... ), , drop=FALSE];
	} else {
		x[ do.call( "order", c( x[ key ], ... ) ), , drop=FALSE ];
	}
}
