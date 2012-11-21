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

## MG_RVG_ALIAS
##
## SUMMARY
##  The Alias method.
##  Given n discrete events with different probabilities P[x] produce a value
##  x consistent with its probability.
##
## DESCRIPTION
## 
##  This algorithmus is based on [1,2]. It is an ingeneous method for
##  generating random variates with finite probability vector which requires
##  a table of size N and needs only one comparision.
##
##  Walker's algorithm does some preprocessing, and provides two array:
##  floating point Q[k] and integer J[k]. A value k is chosen from 0..N-1
##  with equal likelihood, and then a uniform random number u is compared to
##  Q[k].  If it is less than Q[k], then k is returned. Otherwise, J[k] is
##  returned.
##  The method has been generalized in [4] to the "alias-urn" method that
##  works in the exactly same way but uses a larger table of aliases.
##
##  The original algorithm needs 2 uniform random numbers. By reusing only
##  one is necessary (see [6]).
##
##  Walker's original paper describes an O(N^2) algorithm for setting
##  up the Q and J arrays. It had been improved in e.g. [3].
##  This implementation uses Marsaglia's "Robin Hood algorithm" (see [7]):
##  This O(N) algorithm goes through all the p_k's and decides if they are
##  are "poor" or "rich" according to whether they are less than or greater
##  than (>=) the mean value 1 (For convienience we normalize the p_k's,
##  s.t. there sum N instead of 1.). The indices to the poors and the richs
##  are put in separate stacks, and then we work through the stacks together.
##  We take from the next rich on the stack and give to the poor, s.t. it
##  has the average value, and then it is popped from the stack of poors and
##  stored in the tables Q and J.
##  (Walker always wanted to pair up the poorest with the richest.)
##  This reduces the size of the rich and even might become poor, i.e.,
##  it is popped from the stack of richs and pushed on the stack of poors.
##  Since the size of the the two stacks together is <= N, we use one array
##  and store the poors on the beginning and the richs at the and of an
##  array of size N.
##
## SEE ALSO
##  Alastair J. Walker. "An Efficient Method for Generating Discrete Randome Variables with General Distributions", 1977
##  Richard A. Kronmal, Arthur V. Peterson Jr. "On the Alias Method for Generating Random Variable from a Discrete Distribution.", 1979
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

## Generates tables for cutoffs and aliases and returns them in a list object.
mg_rvg_alias.init <- function( p )
{
	qx <- c(); # the cutoffs
	jx <- c(); # the aliases

	rel.tol = .Machine$double.eps^0.5; # square-root of double eps
	n <- length( p ); # len of probs array
	s <- sum( p ); # sum of all probs

	# Creates list of "poor" and "rich" strips
	ratio <- n/s;
	prSize <- n+2;
	poor_rich <- c(); # array of n+2 elements
	poor <- 1;
	rich <- prSize;
	i <- 1;
	#for( i in c(1,n) )
	while ( i <= n )
	{
		qx[i] <- p[i] * ratio; # probability rescaled
		if ( qx[i] >= 1 )
		{
			# rich strip
			poor_rich[ rich ] <- i; # add to list ...
			rich <- rich - 1; # ... and update index
			jx[i] <- i; # init donor (itself)
		} else {
			# poor strip
			poor_rich[ poor ] <- i; # add to list ...
			poor <- poor + 1; # ... and update index
			# (it's not necessary to mark donor)
		}
		i <- i + 1;
	}

	# all other (additional) strips own nothing yet
	while ( i <= n )
	{
		qx[i] <- 0;
		poor_rich[ poor ] <- i;
		poor <- poor + 1;

		i <- i + 1;
	}

	# there must be at least one rich strip
	if ( rich == prSize )
	{
		stop( "[MG::RVG::ALIASINIT] No rich strips found for Robin Hood algorithm." );
	}

	rich <- rich + 1; # must point to the first rich strip yet

	# Makes the "squared histogram" with Robin Hood algorithm (Marsaglia).
	while ( poor != 1 )
	{
		if ( rich > prSize )
		{
			# wrong: assume a neglactable round off error
			break;
		}
		npoor <- poor - 1; # takes the next poor from stack
		jx[ poor_rich[ npoor ] ] <- poor_rich[ rich ]; # store the donor
		qx[ poor_rich[ rich ] ] <- qx[ poor_rich[ rich ] ] - 1 + qx[ poor_rich[ npoor ] ]; # update rich

		# rich might has given too much, so it is poor then
		if ( qx[ poor_rich[ rich ] ] < 1 )
		{
			npoor <- poor_rich[ rich ]; # exchange noveau-poor with former poor in list
			rich <- rich + 1; # remove it from list of rich
		} else {
			poor <- poor - 1; # remove poor from list
		}
	}

	# if there ihas been an round off error, we have to complete the table
	if ( poor != 1 )
	{
		s <- 0; # we estimate the round off error
		while ( poor != 1 )
		{
			npoor <- poor - 1; # take next poor from stack
			s <- s + 1 - qx[ poor_rich[ npoor ] ];
			jx[ poor_rich[ npoor ] ] <- poor_rich[ npoor ]; # mark donor as "not valid"
			qx[ poor_rich[ npoor ] ] = 1; # set probability to 1 (we assume it is very close to one)
			poor <- poor - 1; # remove from list
		}
rel.tol=1e-12
		if ( abs( s ) > rel.tol )
		{
			warning( "[MG::RVG::ALIAS] Sum of deviations very large." );
		}
	}

	return( list( alias = jx, cutoff = qx ) );
}

## Returns a random number for distribution "p" generated with the alias method
## using the given alias-cutoff "ac" parameter.
mg_rvg_alias.rand <- function( n, p, aliases.cutoffs )
{
	if ( missing( aliases.cutoffs ) )
	{
		aliases.cutoffs <- mg_rvg_aliasInit( p );
	}

	r <- c();
	m <- length(p);
	for ( k in 1:n )
	{
		u <- runif(1); # u is U(0,1)
		u <- u * m; # u is U(0,n)
		i <- ceiling( u ); # i is U({1,2,...,n})
		u <- i - u; # reuse rand number

		if ( u <= aliases.cutoffs$cutoff[i] )
		{
			r[k] <- i;
		} else {
			r[k] <- aliases.cutoffs$alias[i];
		}
	}

	return( r );
}
