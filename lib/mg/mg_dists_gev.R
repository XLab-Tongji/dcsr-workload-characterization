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

## MG_DISTS_GEV.R
##
## SUMMARY
##  Generalized Extreme Value distribution.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

mg_dists_dgev <- function( x, loc = 0, scale = 1, shape = 0, log = FALSE )
{
	haveevd <- require( "evd" );

	if ( haveevd )
	{
		return( evd::dgev( x, loc, scale, shape, log ) );
	}

	stop( "Not yet implemented!" ); #FIXME 
}

mg_dists_pgev <- function( q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE )
{
	haveevd <- require( "evd" );

	if ( haveevd )
	{
		p <- evd::pgev( q, loc, scale, shape, lower.tail );
		if ( log.p )
		{
			p <- log( p );
		}
		return( p );
	}

	stop( "Not yet implemented!" ); #FIXME 
}

mg_dists_qgev <- function( p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE )
{
	haveevd <- require( "evd" );

	if ( haveevd )
	{
		if ( log.p )
		{
			p <- exp( p );
		}
		return( evd::qgev( p, loc, scale, shape, lower.tail ) );
	}

	stop( "Not yet implemented!" ); #FIXME 
}

mg_dists_rgev <- function( n, loc = 0, scale = 1, shape = 0 )
{
	haveevd <- require( "evd" );

	if ( haveevd )
	{
		return( evd::rgev( n, loc, scale, shape ) );
	}

	stop( "Not yet implemented!" ); #FIXME 
}
