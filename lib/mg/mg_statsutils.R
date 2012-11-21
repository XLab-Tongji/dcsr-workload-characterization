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

## MG_STATSUTILS
##
## SUMMARY
##  Some useful statistics utilities
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

## MG_STATSUTILS_ZSCORE(x,m=mean(x),s=sd(x))
##
## SUMMARY
##  Returns the z-score value.
##
## DESCRIPTION
##  Standardize the given value respect to the Standard Normal distribution
##   N(0,1).
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##
mg_statsutils_zscore <- function( x, m = mean(x), s = sd(x) )
{
	return( ( x - m ) / s );
}

## MG_STATSUTILS_BOOTSTRAP(x,nboot,theta)
##
## SUMMARY
##  Generates @c nboot resamples and apply @c theta estimator
##
## E.g.: generate 1000 bootstrap sample medians:
##  > data <- c( 3.19,4.26,4.47,4.53,4.67,4.69,5.78,6.79,9.37,12.75 )
##  > bootsamples <- mg_statsutils_bootstrap( data, 1000, "median")
##
mg_statsutils_bootstrap <- function( x, nboot, theta )
{
	data <- matrix( sample( x, size = length(x)*nboot, replace = TRUE), nrow = nboot );
	return( apply( data, 1, theta ) );
}

## MG_STATSUTILS_CV
##
## SUMMARY
##  Coefficient of variation.
##
## PARAMS
##  x: data values.
##
## RETURN
##  The coefficient of variation.
##
## DESCRIPTION
##  Computes the coefficient of variation statistic, defined as CV = sigma / mu,
##  where 'sigma' is the standard error and 'mu' is the mean.
##
mg_statsutils_cv <- function( x )
{
	return( sd( x ) / mean( x ) );
}

## MG_STATSUTILS_PIT
##
## SUMMARY
##  Probability Intetral Transformation (aka Uniform Transformation).
##
mg_statsutils_pit <- function( x, dist, ... )
{
	return( mg_dists_cdf( dist, x, ... ) );
}

## MG_STATSUTILS_CONFINT.MEAN(x,level=0.95,side=c("two.sided", "less", "greater"), s = NULL )
##
## SUMMARY
##  Confidence Interval for mean
##
mg_statsutils_confint.mean <- function( x, level = 0.95, side = c("two.sided", "less", "greater"), s = NULL )

{
	side <- match.arg( side );

	s.unknown <- FALSE;
	if ( is.null(s) )
	{
		s <- sd(x);
		s.uknown <- TRUE;
	}

	x.bar <- mean(x);
	a <- (1 - level);
	n <- length(x);

	l <- -Inf;
	u <- +Inf;

	if ( side == "two.sided" )
	{
		z <- NaN;
		if ( s.unknown )
		{
			z <- qt( a/2, df = n-1, lower.tail = FALSE );
		} else {
			z <- qnorm( a/2, lower.tail = FALSE );
		}
		zz <- z * s/sqrt(n);
		l <- x.bar - zz;
		u <- x.bar + zz;
	} else {
		z <- NaN;
		if ( s.unknown )
		{
			z <- qt( a, df = n-1 );
		} else {
			z <- qnorm( a );
		}
		zz <- z * s/sqrt(n);
		if ( side == "less" )
		{
			l <- x.bar - zz;
		} else {
			u <- x.bar + zz;
		}
	}

	r <- as.data.frame( list( mean = x.bar, lower = l, upper = u ) );
	class(r) <- "table";

	return( r );
}

## MG_STATSUTILS_CONFINT.VAR(x,level=0.95,side=c("two.sided", "less", "greater"))
##
## SUMMARY
##  Confidence Interval for variance.
##
mg_statsutils_confint.var <- function( x, level = 0.95, side = c("two.sided", "less", "greater") )

{
	side <- match.arg( side )

	s2 <- var(x);

	a <- (1 - level);
	n <- length(x);

	l <- -Inf;
	u <- +Inf;

	ss2 <- (n-1) * s2;

	if ( side == "two.sided" )
	{
		l <- ss2 / qchisq( a/2, df = n-1, lower.tail = FALSE );
		u <- ss2 / qchisq( 1-a/2, df = n-1, lower.tail = FALSE );
	} else {
		if ( side == "less" )
		{
			l <- ss2 / qchisq( a, df = n-1, lower.tail = FALSE );
		} else {
			u <- ss2 / qchisq( 1-a, df = n-1, lower.tail = FALSE );
		}
	}

	r <- as.data.frame( list( var = s2, lower = l, upper = u ) );
	class(r) <- "table";

	return( r );
}

## MG_STATSUTILS_CONFINT.MEDIAN(x,level=0.95,side=c("two.sided", "less", "greater"))
##
## SUMMARY
##  Confidence Interval for median.
##
## SEE ALSO
##  Lehmann (Nonparametrics: Statistical Methods Based on Ranks, Holden-Day, 1975, p. 182-183).
##  (for alternatives, see http://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/mediancl.htm)
##
mg_statsutils_confint.median <- function( x, level = 0.95, side = c("two.sided", "less", "greater") )
{
	## >>> TODO: add one-sided (less, greater) intervals <<<

	side <- match.arg( side )

	a <- 1 - level;

	v <- sort(x, na.last = NA);
	n <- length(x);
	if( n > 0 )
	{
		m <- median(x);

		if ( side == "two.sided" )
		{
			#i <- qbinom(0.025, n, 0.5);
			i <- qbinom(a/2, n, a);
			if(i > 0)
			{
				r <- c(m, v[i], v[n - i + 1]);
			} else {
				r <- c(m, NA, NA);
			}
		} else {
			r <- c(m, NA, NA);
			warning( "[MG::STATSUTILS::CONFINT::MEDIAN] One-sided confidence intervals not yet implemented." );
		}
	} else {
		r <- c(NA, NA, NA);
	}

	r <- as.data.frame( list( median = r[1], lower = r[2], upper = r[3] ) );
	class(r) <- "table";

	return( r );
}

## MG_STATSUTILS_PERCENTILE
##
## SUMMARY
##  Gives the percentiles of the sample in X.
##  Y = PRCTILE(X,P) returns a value that is greater than P percent
##  of the values in X. For example, if P = 50  Y is the median of X. 
##
##  P may be either a scalar or a vector. For scalar P, Y is a row   
##  vector containing Pth percentile of each column of X. For vector P,
##  the ith row of Y is the P(i) percentile of each column of X.
##
#mg_statsutils_percentile <- function( x, p )
#{
#
#	if ( any(p > 100) | any(p < 0) )
#	{
#		stop( "P must take values between 0 and 100" );
#	}
#
#	xx = sort(x);
#	xx.len = length(xx);
#
#	if ( xx.len == 1 )
#	{
#		q = 100*(0.5:x.len - 0.5)./x.len;
#		xlims <- range(xx);
#		xx = [min(xx); xx(:); max(xx)];
#	} else {
#		q = 100*(0.5:m - 0.5)./m;
#		xx = [min(x); xx; max(x)];
#	}
#
#	q = [0 q 100];
#
#	y = interp1(q,xx,p);
#}
