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

##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

## Z-statistic
##   Given an n-sample x1, x2, ..., xn, and the "real" standard deviation (i.e.
##   standard deviation coming from the underlying population), let
##     sm = (1/n) " sum(i=1 to n; xi) [sample mean]
##   then:
##     z = ( sm - m ) / ( s / sqrt( n ) )
mg_zstat <- function( x, m, s )
{
	z <- sqrt( length( x ) ) * ( mean( x ) - m ) / s;

	# return value
        z;
}

## Two-sided z-test
##   H0: mu = mu0
##   H1: mu <> mu0
##
## Input Parameters:
##   x: a sample.
##   s: the "real" (i.e. population) standard deviation.
##   mu0: the parameter value used for the two-side test (i.e. H0: mu = mu0).
##   mu1: a parameter value for the alternative hypothesis (i.e. H1: mu = mu1).
##   a: significance level.
##
## Output Values:
##   rejected: TRUE if H0 has been rejected.
##   pval: p-value; is the smallest significance level for which H0 is rejected,
##         i.e. H0 is rejected for all signficance levels a such that a >= pval.
##   a: the type I error; it is the probability of not rejecting H0 when it is
##      true (is equal to the significance level).
##   b: the type II error; it is the probability of not rejecting H0 when it is
##      false (i.e. H1 is true).
##   xleft: left critical value, i.e. the left limit for a (1-a) confidence
##          interval.
##   xright: right critical value, i.e. the right limit for a (1-a) confidence
##           interval.
##   power: the power of the test; this is the probability of rejecting H0 if H1
##          is true.
mg_ztest2side <- function( x, s, mu0, mu1, a = 0.05 )
{
	# Initializations
	n <- length( x );
	sn <- s / sqrt( n ); # standard deviation for sample mean

	# Output vals
	rejected <- pval <- b <- xcleft <- xcright <- b <- NULL;

	z0 <- mg_zstat( x, mu0, s ); # Z-statistic
	pval <- 2 * pnorm( -abs( z0 ) ); # p-value
	#if ( abs( z0 ) > qnorm( 1-a/2) ) # Alternative #1
	#if ( abs( z0 ) > qnorm( a/2, lower.tail=FALSE) ) # Alternative #2
	if ( a >= pval ) # Alternative #3
	{
		# H0 will be rejected for all significance levels >= pval
		reject <- TRUE;
	} else {
		reject <- FALSE;
	}

	# Calculates Critical values and Type II error
	err <- qnorm( 1 - a/2 ) * sn;
	xcleft <- mu0 - err; # Left critical value
	xcright <- mu0 + err; # Right critical value
	if ( !missing( mu1 ) )
	{
		#b <- pnorm( (mu0-mu1)/sn + qnorm(1-a/2) ) - pnorm( (mu0-mu1)/sn - qnorm(1-a/2) ); # Type II error
		zcleft <- mg_zstat( xcleft, mu1, sn );
		zcright <- mg_zstat( xcright, mu1, sn );
		b <- pnorm( zcright ) - pnorm( zcleft ); # Type II error
	}

	# return values
	list(
		z = z0, # z-statistic
		rejected = reject, # TRUE if null hypothesis is rejected
		pval = pval, # p-value
		a = a, # Type I error
		b = b, # Type II error
		xcleft = xcleft, # Left critical value
		xcright = xcright, # Right critical value
		power = ifelse( !is.null( b ), 1 - b, 0 ) # Power of the test
	);
}

## One-sided (less than) z-test
##   H0: mu >= mu0 (or mu = mu0)
##   H1: mu < mu0
mg_ztest1sidelt <- function( x, s, mu0, mu1, a = 0.05 )
{
	# Initializations
	n <- length( x );
	sn <- s / sqrt( n ); # standard deviation for sample mean

	# Output vals
	rejected <- pval <- b <- xc <- b <- NULL;

	z0 <- mg_zstat( x, mu0, s ); # Z-statistic
	pval <- pnorm( z0 ); # p-value [ P(Z <= z0) ]
	#if ( z0 < qnorm( 1-a ) ) # Alternative #1
	#if ( z0 < qnorm( a, lower.tail=FALSE) ) # Alternative #2
	if ( a >= pval ) # Alternative #3
	{
		# H0 will be rejected for all significance levels >= pval
		reject <- TRUE;
	} else {
		reject <- FALSE;
	}

	# Calculates Critical values and Type II error
	err <- qnorm( 1 - a ) * sn;
	xc <- mu0 - err; # Critical value
	if ( !missing( mu1 ) )
	{
		#b <- pnorm( (mu0-mu1)/sn - qnorm(1-a/2) ); # Type II error
		zc <- mg_zstat( xc, mu1, sn );
		b <- pnorm( zc ); # Type II error
	}

	# return values
	list(
		z = z0, # z-statistic
		rejected = reject, # TRUE if null hypothesis is rejected
		pval = pval, # p-value
		a = a, # Type I error
		b = b, # Type II error
		xc = xc, # Critical value
		power = ifelse( !is.null( b ), 1 - b, 0 ) # Power of the test
	);
}

## One-sided (greater than) z-test
##   H0: mu <= mu0 (or mu = mu0)
##   H1: mu > mu0
mg_ztest1sidegt <- function( x, s, mu0, mu1, a = 0.05 )
{
	# Initializations
	n <- length( x );
	sn <- s / sqrt( n ); # standard deviation for sample mean

	# Output vals
	rejected <- pval <- b <- xc <- b <- NULL;

	z0 <- mg_zstat( x, mu0, s ); # Z-statistic
	#pval <- pnorm( -z0 ); # p-value [ P( Z <= -z0 ) = P(Z >= z0) ]
	pval <- pnorm( z0, lower.tail = FALSE ); # p-value [ P(Z >= z0) ]
	#if ( z0 > qnorm( 1-a ) ) # Alternative #1
	#if ( z0 > qnorm( a, lower.tail=FALSE) ) # Alternative #2
	if ( a >= pval ) # Alternative #3
	{
		# H0 will be rejected for all significance levels >= pval
		reject <- TRUE;
	} else {
		reject <- FALSE;
	}

	# Calculates Critical values and Type II error
	err <- qnorm( 1 - a ) * sn;
	xc <- mu0 + err; # Critical value
	if ( !missing( mu1 ) )
	{
		#b <- pnorm( (mu0-mu1)/sn + qnorm(1-a/2) ); # Type II error
		zc <- mg_zstat( xc, mu1, sn );
		b <- pnorm( zc ); # Type II error
	}

	# return values
	list(
		z = z0, # z-statistic
		rejected = reject, # TRUE if null hypothesis is rejected
		pval = pval, # p-value
		a = a, # Type I error
		b = b, # Type II error
		xc = xc, # Critical value
		power = ifelse( !is.null( b ), 1 - b, 0 ) # Power of the test
	);
}
