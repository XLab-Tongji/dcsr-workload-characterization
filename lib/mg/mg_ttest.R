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

## T-statistic
##   Given an n-sample x1, x2, ..., xn, let
##     sm = (1/n) " sum(i=1 to n; xi) [sample mean]
##   then:
##     t = ( sm - m ) / ( s / sqrt( n ) )
mg_tstat <- function( x, m, s )
{
	if ( missing( s ) )
	{
		s <- sd( x );
	}

	t <- ( mean( x ) - m ) / ( s / sqrt( length( x ) ) );

	# return value
	t;
}

## Two-sided t-test
##   H0: mu = mu0
##   H1: mu <> mu0
##
## Input Parameters:
##   x: a sample.
##   mu0: the parameter value used for the two-side test (i.e. H0: mu = mu0).
##   mu1: a parameter value for the alternative hypothesis (i.e. H1: mu = mu1).
##   a: significance level.
##
## Return Value:
##   A "list" with following fields:
##   * rejected: TRUE if H0 has been rejected.
##   * pval: p-value; is the smallest significance level for which H0 is
##           rejected, i.e. H0 is rejected for all signficance levels a such
##           that a >= pval.
##   * a: the type I error; it is the probability of not rejecting H0 when it is
##        true (is equal to the significance level).
##   * b: the type II error; it is the probability of not rejecting H0 when it
##        is false (i.e. H1 is true).
##   * xleft: left critical value, i.e. the left limit for a (1-a) confidence
##            interval.
##   * xright: right critical value, i.e. the right limit for a (1-a) confidence
##             interval.
##   * power: the power of the test; this is the probability of rejecting H0 if
##            H1 is true.
mg_ttest2side <- function( x, mu0, mu1, s, a = 0.05)
{
	# Initializations
	n <- length( x );
	df <- n - 1; # degree of freedom
	if ( missing( s ) )
	{
		s <- sd( x ); # sample standard deviation
	}
	sn <- s / sqrt( n ); # standard deviation estimator of sample mean

	# Output vals
	rejected <- pval <- b <- xcleft <- xcright <- b <- NULL;

	t0 <- mg_tstat( x, mu0, s ); # T-statistic
	#pval <- 2 * pt( abs( t0 ), df, lower.tail = FALSE ); # p-value
	pval <- 2 * pt( -abs( t0 ), df ); # p-value
	#if ( abs( t0 ) > qt( 1 - a/2, df ) ) # Alternative #1
	#if ( abs( t0 ) > qt( 1 - a/2, df ) ) # Alternative #2
	if ( a >= pval ) # Alternative #3
	{
		# H0 will be rejected for all significance levels >= pval

		reject <- TRUE;
	} else {
		reject <- FALSE;
	}

	# Calculates Critical values and Type II error
	err <- qt( 1 - a/2, df ) * sn;
	xcleft <- mu0 - err; # Left critical value
	xcright <- mu0 + err; # Right critical value
	if ( !missing( mu1 ) )
	{
		#b <- pt( (mu0-mu1)/sn + qt(1-a/2, df), df ) - pt( (mu0-mu1)/sn - qt(1-a/2, df), df ); # Type II error
		tcleft <- mg_tstat( xcleft, mu1, sn );
		tcright <- mg_tstat( xcright, mu1, sn );
		b <- pt( tcright, df ) - pt( tcleft, df ); # Type II error
	}

	# return values
	list(
		t = t0, # t-statistic
		rejected = reject, # TRUE if null hypothesis is rejected
		pval = pval, # p-value
		a = a, # Type I error
		b = b, # Type II error
		xcleft = xcleft, # Left critical value
		xcright = xcright, # Right critical value
		power = ifelse( !is.null( b ), 1 - b, 0 ) # Power of the test
	);
}

## One-sided (less than) t-test.
##   H0: mu = mu0 (or mu = mu0)
##   H1: mu < mu0
mg_ttest1sidelt <- function( x, mu0, mu1, s, a = 0.05 )
{
	# Initializations
	n <- length( x );
	df <- n - 1; # degree of freedom
	if ( missing( s ) )
	{
		s <- sd( x ); # sample standard deviation
	}
	sn <- s / sqrt( n ); # standard deviation estimator of sample mean

	# Output vals
	rejected <- pval <- b <- xc <- b <- NULL;

	t0 <- mg_tstat( x, mu0, s ); # T-statistic
	#pval <- pt( -t0, df, lower.tail = FALSE ); # p-value
	pval <- pt( t0, df ); # p-value [ P(T <= t0 )
	#if ( t0 < -qt( 1 - a, df ) ) # Alternative #1
	#if ( t0 < -qt( 1 - a, df ) ) # Alternative #2
	if ( a >= pval ) # Alternative #3
	{
		# H0 will be rejected for all significance levels >= pval

		reject <- TRUE;
	} else {
		reject <- FALSE;
	}

	# Calculates Critical values and Type II error
	err <- qt( 1 - a, df ) * sn;
	xc <- mu0 + err; # Critical value
	if ( !missing( mu1 ) )
	{
		#b <- pt( (mu0-mu1)/sn + qt(1-a/2, df), df ); # Type II error
		tc <- mg_tstat( xc, mu1, sn );
		b <- pt( tc, df ); # Type II error
	}

	# return values
	list(
		t = t0, # t-statistic
		rejected = reject, # TRUE if null hypothesis is rejected
		pval = pval, # p-value
		a = a, # Type I error
		b = b, # Type II error
		xc = xc, # Critical value
		power = ifelse( !is.null( b ), 1 - b, 0 ) # Power of the test
	);
}

## One-sided (greater than) t-test.
##   H0: mu <= mu0 (or mu = mu0)
##   H1: mu > mu0
mg_ttest1sidegt <- function( x, mu0, mu1, s, a = 0.05 )
{
	# Initializations
	n <- length( x );
	df <- n - 1; # degree of freedom
	if ( missing( s ) )
	{
		s <- sd( x ); # sample standard deviation
	}
	sn <- s / sqrt( n ); # standard deviation estimator of sample mean

	# Output vals
	rejected <- pval <- b <- xc <- b <- NULL;

	t0 <- mg_tstat( x, mu0, s ); # T-statistic
	pval <- pt( t0, df, lower.tail = FALSE ); # p-value
	#pval <- pt( -t0, df ); # p-value
	#if ( t0 > qt( 1 - a, df ) ) # Alternative #1
	#if ( t0 > qt( 1 - a, df ) ) # Alternative #2
	if ( a >= pval ) # Alternative #3
	{
		# H0 will be rejected for all significance levels >= pval

		reject <- TRUE;
	} else {
		reject <- FALSE;
	}

	# Calculates Critical values and Type II error
	err <- qt( 1 - a, df ) * sn;
	xc <- mu0 - err; # Critical value
	if ( !missing( mu1 ) )
	{
		#b <- pt( (mu0-mu1)/sn - qt(1-a/2, df), df ); # Type II error
		tc <- mg_tstat( xc, mu1, sn );
		b <- pt( tc, df ); # Type II error
	}

	# return values
	list(
		t = t0, # t-statistic
		rejected = reject, # TRUE if null hypothesis is rejected
		pval = pval, # p-value
		a = a, # Type I error
		b = b, # Type II error
		xc = xc, # Critical value
		power = ifelse( !is.null( b ), 1 - b, 0 ) # Power of the test
	);
}
