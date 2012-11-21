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

## MG_FIT
##
## SUMMARY
##  A collection of fitting methods for several distribution.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

source( "lib/mg/mg_consts.R" );
source( "lib/mg/mg_dists.R" );
source( "lib/mg/mg_moments.R" );
source( "lib/mg/mg_fit_utils.R" );

mg_fit_mleBeta <- function( x )
{
	library( MASS );

	if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
	{
		stop("[MG::FIT::MLEBETA] 'x' must be a non-empty numeric vector");
	}

	# BEGIN normaliza x in the range [0,1]
	x <- x[ x > 0 & is.finite(x) ]; # remove all zeros
	if ( any( x > 1 ) )
	{
		# normalize between [0,1]
		x <- x / sum(x);
	}
	# END normaliza x in the range [0,1]

	# BEGIN initial values calculation taken from MATLAB "betafit" function
	n <- length( x );
	tmp1 <- prod( (1-x)^(1/n) );
	tmp2 <- prod( x^(1/n) );
	tmp3 <- 1 - tmp1 - tmp2;
	s1 <- 0.5*(1-tmp1) / tmp3;
	s2 <- 0.5*(1-tmp2) / tmp3;
	# END initial values calculation taken from MATLAB "betafit" function

	x <- MASS::fitdistr( x, "beta", list( shape1 = s1, shape2 = s2 ) );

	return(
		list(
			shape1 = as.numeric( x$estimate["shape1"] ),
			shape2 = as.numeric( x$estimate["shape2"] )
		)
	);
}

## Estimates Beta parameters (shape1 and shape2) with the method-of-moments.
mg_fit_momBeta <- function( x )
{
	# BEGIN normaliza x in the range [0,1]
	x <- x[ x > 0 & is.finite(x) ]; # remove all zeros
	if ( any( x > 1 ) )
	{
		# normalize between [0,1]
		x <- x / sum(x);
	}
	# END normaliza x in the range [0,1]

	m <- mean( x );
	v <- var( x );
	s1 <- m * ( m * (1-m)/v - 1);
	s2 <- (1-m) * ( m * (1-m)/v - 1);

	return( list( shape1 = s1, shape2 = s2 ) );
}

mg_fit_mleBinomial <- function( x )
{
	library( MASS );

	if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
	{
		stop("[MG::FIT::MLEBINOMIAL] 'x' must be a non-empty numeric vector");
	}

	# Ensures x is in the set {0,...,size}
	x <- x[ x >= 0 & is.finite(x) ];

	x <- MASS::fitdistr( x, "binomial" );

	return(
		list(
			size = as.numeric( x$estimate["size"] ),
			prob = as.numeric( x$estimate["prob"] )
		)
	);
}

mg_fit_mleCauchy <- function( x )
{
	library( MASS );

	if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
	{
		stop("[MG::FIT::MLECAUCHY] 'x' must be a non-empty numeric vector");
	}

	x <- MASS::fitdistr( x, "cauchy" );

	return(
		list(
			location = as.numeric( x$estimate["location"] ),
			scale = as.numeric( x$estimate["scale"] )
		)
	);
}

mg_fit_mleChiSquared <- function( x )
{
	library( MASS );

	if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
	{
		stop("[MG::FIT::MLECHISQUARED] 'x' must be a non-empty numeric vector");
	}

	# Ensures x is in the range (0,+Inf)
	x <- x[ x > 0 & is.finite(x) ];
	x <- MASS::fitdistr( x, "chi-squared", start = list( df = 1 ) );

	return(
		list(
			df = as.numeric( x$estimate["df"] )
		)
	);
}

## Fit a distribution with a Phase-Type using the moment matching method.
mg_fit_momCPH <- function( x )
{
	mu1 <- mean( x );
	mu2 <- mg_moments_moment( x, 2, central = FALSE );
	mu3 <- mg_moments_moment( x, 3, central = FALSE );
#TODO

	#cat( paste( "octave -q --eval 'addpath(\"", getwd(), "/../MATLAB\"); [tau,T] = matching3PH([", formatC(mu1, format = "f"), ",", formatC(mu2, format = "f"), ",", formatC(mu3, format = "f"), "])'", sep = "" ), "\n" );
#warning( paste("PHASE-TYPE FIT: MU1: ", mu1, " - MU2: ", mu2, " - MU3: ", mu3, sep = "" ), immediate = TRUE );#XXX
	res <- try( system( paste( "octave -q --eval 'addpath(\"", getwd(), "/../MATLAB\"); [tau,T] = matching3PH([", formatC(mu1, format = "f"), ",", formatC(mu2, format = "f"), ",", formatC(mu3, format = "f"), "]); format free; disp([tau;T])'", sep = "" ), TRUE ), TRUE );
#print( "PHASE-TYPE FIT: RES: " );#XXX
#print( res );#XXX
	rows <- unlist( strsplit( res, "\n" ) );
	rows <- rows[ nchar( rows ) > 0 ]; # removes blank rows
#print( "PHASE-TYPE FIT: ROWS: " );#XXX
#print( rows, immediate = TRUE );#XXX
	nums <- unlist( strsplit( rows[1], "[[:space:]]" ) );
	nums <- as.numeric( nums[ nchar( nums ) > 0 ] );
#print( "PHASE-TYPE FIT: NUMS: " );#XXX
#print( nums );#XXX
	alpha <- as.numeric( nums[ nchar( nums ) > 0 ] );
#print( paste("PHASE-TYPE FIT: ALPHA: ", length(alpha), sep = "" ) );#XXX
	#nums <- unlist( strsplit( res, "[[:space:]]" ) );
	nums <- unlist( strsplit( rows[-1], "[[:space:]]" ) );
	nums <- as.numeric( nums[ nchar( nums ) > 0 ] );
#print( "PHASE-TYPE FIT: NUMS: " );#XXX
#print( nums );#XXX
	T <- matrix( nums, nrow = length(rows)-1, ncol = length(rows)-1, byrow = TRUE );
#	n <- length( res );
#warning( paste("PHASE-TYPE FIT: N: ", n, sep = "" ), immediate = TRUE );#XXX
	##alpha <- nums[ 1:(n-1) ];
	#alpha <- nums[ 1:(n-1) ];
	#T <- matrix( nums[ -(1:(n-1)) ], nrow = (n-1), ncol = (n-1), byrow = TRUE );
	#T <- matrix( nums[ -(1:(n-2)) ], nrow = (n-2), ncol = (n-2), byrow = TRUE );
#warning( paste("PHASE-TYPE FIT: T: ", dim(T), sep = "" ), immediate = TRUE );#XXX
	# BEGIN Overall interarrival time
	#alpha = c( 1, 0 );
	#T = matrix( c( -0.49636, 0.22716, 0.00000, -0.15057 ), nrow = 2, ncol = 2, byrow = TRUE );
	# END Overall interarrival time

	return( list( alpha = alpha, T = T ) );
}

mg_fit_mleExponential <- function( x )
{
	library( MASS );

	#dots <- names(list(...))
	#dots <- dots[!is.element(dots, c("upper", "lower"))]
	if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
	{
		stop("[MG::FIT::MLEEXPONENTIAL] 'x' must be a non-empty numeric vector");
	}

	# Ensures x is in the range [0,+Inf)
	x <- x[ x >= 0 & is.finite(x) ];
	x <- MASS::fitdistr( x, "exponential" );

	return(
		list(
			rate = as.numeric( x$estimate["rate"] )
		)
	);
}

mg_fit_mleFisherF <- function( x )
{
	library( MASS );

	if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
	{
		stop("[MG::FIT::MLEFISHERF] 'x' must be a non-empty numeric vector");
	}

	# Ensures x is in the range (0,+Inf)
	x <- x[ x > 0 & is.finite(x) ];
	x <- MASS::fitdistr( x, "f", start = list( df1 = 1, df2 = 1 ) );

	return(
		list(
			df1 = as.numeric( x$estimate["df1"] ),
			df2 = as.numeric( x$estimate["df2"] )
		)
	);
}

#mg_fit_mleFrechet <- function( x )
#{
#	return( mg_dists_frechet.fit.mle( x ) );
#}

mg_fit_mleGamma <- function( x )
{
	library( MASS );

	if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
	{
		stop("[MG::FIT::MLEGAMMA] 'x' must be a non-empty numeric vector");
	}

# FIXME: fitdistr throws error (even with x[x>0])
#	#x <- x[ x > 0 ];
#	x <- MASS::fitdistr( x, "gamma" );
#
#	return( list( shape = x$estimate["shape"], rate = x$estimate["rate"], scale = x$estimate["scale"] ) );

	# The Gamma distribution has support in [0,+Inf)
	x <- x[ x >= 0 & is.finite(x) ];

	if ( any( x == 0 ) )
	{
		# MLE not possible => use method of moments
		warning( "[MG::FIT::MLEGAMMA] Zeros in data -- using method of moments as estimators." );
		return( mg_fit_momGamma( x ) );

	}

	fit <- try( MASS::fitdistr( x, "gamma" ) );
	if ( is( fit, "try-error" ) )
	{
		# Use an approximation based on Newton method for shape parameter
		warning( "[MG::FIT::MLEGAMMA] MLE failed -- using Newton approximation." );

		n <- length( x );
		s <- log(mean(x)) - mean(log(x)); 
		shape <- (3-s+sqrt((s-3)^2 + 24*s))/(12*s);
		scale <- mean( x ) / shape;

		fit <- list( shape = shape, scale = scale );
	} else {
		fit <- list(
			shape = as.numeric( fit$estimate["shape"] ),
			scale = as.numeric( fit$estimate["scale"] )
		);
	}

	return( fit );
}

mg_fit_momGamma <- function( x )
{
	# Ensures x is in the range [0,+Inf)
	x <- x[ x >= 0 & is.finite(x) ];

	xbar <- mean(x);
	s2 <- var(x);
	shapehat <- xbar * xbar / s2;
	ratehat <- shapehat / xbar; # i.e. xbar / s2

	return( list( shape = shapehat, rate = ratehat ) );
}

mg_fit_mleGEV <- function( x )
{
	library( evd );

	fit <- evd::fgev( x );

	return(
		list(
			loc = as.numeric( fit$estimate["loc"] ),
			scale = as.numeric( fit$estimate["scale"] ),
			shape = as.numeric( fit$estimate["shape"] )
		)
	);
}

##mg_fit_mleGPD <- function( x, loc = min(x) )
#mg_fit_mleGPD <- function( x, loc = 0 )
#{
#	return( mg_dists_gpd.fit.mle( x, loc = loc ) );
#}

#mg_fit_mleGumbel <- function( x )
#{
#	return( mg_dists_gumbel.fit.mle( x ) );
#}

mg_fit_emHyperExp <- function( x )
{
	#FIXME!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	return( list( probs = c( 0.6966325, 0.3033675 ), rates = c( 0.1508341, 0.6963743 ) ) ); # LCG-2005-0 Overall Interarrival times
	#return( list( probs = c( 0.3014699, 0.6985301 ), rates = c( 0.07248913, 0.07248913 ) ) ); # LCG-2005-0 VO lhcb Interarrival times
}

mg_fit_mleLogistic <- function( x )
{
	library( MASS );

	if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
	{
		stop("[MG::FIT::MLELOGISTIC] 'x' must be a non-empty numeric vector");
	}

	x <- MASS::fitdistr( na.omit(x), "logistic" );

	return(
		list(
			location = as.numeric( x$estimate["location"] ),
			scale = as.numeric( x$estimate["scale"] )
		)
	);
}

#FIXME: improve handling zero or negative values
mg_fit_mleLogNormal <- function( x )
{
	library( MASS );

	if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
	{
		stop("[MG::FIT::MLELOGNORMAL] 'x' must be a non-empty numeric vector");
	}

#	# sanity checks
#	if ( any( oops <- ( x <= 0) ) )
#	{
#		#x[ oops ] <- -Inf;
#		x[ oops ] <- NA;
#	}

	# Ensures x is in the range (0,+Inf)
	x <- x[ x > 0 & is.finite(x) ];
	x <- MASS::fitdistr( na.omit(x), "log-normal" );

	return(
		list(
			meanlog = as.numeric( x$estimate["meanlog"] ),
			sdlog = as.numeric( x$estimate["sdlog"] )
		)
	);
}

mg_fit_emMMPP <- function( x )
{
	stop( "[MG::FIT::EMMMPP] Not yet implemented" );
}

mg_fit_mleNormal <- function( x )
{
	library( MASS );

	#dots <- names(list(...))
	#dots <- dots[!is.element(dots, c("upper", "lower"))]
	if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
	{
		stop("[MG::FIT::MLENORMAL] 'x' must be a non-empty numeric vector");
	}

	x <- MASS::fitdistr( x, "normal" );

	return(
		list(
			mean = as.numeric( x$estimate["mean"] ),
			sd = as.numeric( x$estimate["sd"] )
		)
	);
}

##
## Log-Likelihood function L(a, xmin) for Pareto:
##  L(a, xmin) = prod{ i=1, n, k * (xmin^a / x_i^(a+1) } = a^n * min^(na) * prod{ i=1, n, 1 / x_i^(k+1) } 
## (see http://en.wikipedia.org/wiki/Pareto_distribution)
mg_fit_mlePareto <- function( x )
{
	#params <- mg_fit_hillEstimator( x );
	#params <- mg_fit_yuleEstimator( x, xmin = xmin );
#print( params );#XXX
	#return ( params );

	x <- x[ is.finite(x) ];

	# Estimates the Xmin (FIXME: maybe to be enhanced)
	xm.est <- max( min( x ), .Machine$double.eps^.25 ); # this is the estimation of xmin (avoid xmin==0)

	# since we use logarithms (log(x)) we throw away values <= 0
	x <- x[ x > 0 ];

	# Starting value for Xmin; note we don't use xm.est since it can be = 0
	xm <- min( x );
	n <- length( x );
	alpha <- n / ( sum( log( x ) ) - n * log( xm ) );

	return( list( shape = alpha,  xmin = xm.est ) );
}

mg_fit_hillPareto <- function( x )
{
	# Ensure x is in safe range
	x <- x[ x > 0 & is.finite(x) ];

	# this is the estimation of xmin (avoid xmin==0)
	#xm <- min( x );
	xm <- max( min( x ), .Machine$double.eps^.25 );

	hill <- mg_fit_hillEstimator( x, length( x ) );
	alpha <- hill$statistic[ length( hill$statistic ) ]; # get the last hill estimate
	alpha <- 1 / alpha; # the hill estimator is the inverse tail index

	return( list( shape = alpha,  xmin = xm ) );
}

mg_fit_yulePareto <- function( x )
{
	# Ensure x is in safe range
	x <- x[ x > 0 & is.finite(x) ];

	# this is the estimation of xmin (avoid xmin==0)
	xm <- min( x );
	xm <- max( min( x ), .Machine$double.eps^.25 );

	res <- mg_fit_yuleEstimator( x, xmin = xm );

	return( list( shape = res$alpha,  xmin = res$xmin ) );
}

mg_fit_mlePoisson <- function( x )
{
	library( MASS );

	#dots <- names(list(...))
	#dots <- dots[!is.element(dots, c("upper", "lower"))]
	if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
	{
		stop("[MG::FIT::MLEPOISSON] 'x' must be a non-empty numeric vector");
	}

	# Ensures x is in the set {0,1,...}
	x <- x[ x >= 0 & is.finite(x) ];

	x <- MASS::fitdistr( x, "Poisson" );

	return(
		list(
			lambda = as.numeric( x$estimate["lambda"] )
		)
	);
}

mg_fit_mleReverseWeibull <- function( x )
{
	stop( "[MG::FIT::MLEREVERSEWEIBULL] Not yet implemented" );
}

## Fits Skewed Stable distributions with MLE method.
mg_fit_mleSkewedStable <- function( x )
{
	return( mg_fit_mleStable( x ) );
}

## Fits Symmetric Stable distributions with MLE method.
mg_fit_mleSymStable <- function( x )
{
	res <- mg_fit_mleStable( x );

	# removes unused information
	res$beta <- res$gamma <- res$delta <- NULL;

	return( res );
}

## Fits Stable distributions with MLE method.
#mg_fit_mleStable <- function( x )
#{
#	return( mg_dists_stable.fit.mle( x, doplot = FALSE ) );
#}

## Fits Skewed Stable distributions with McCulloch's Quantile method.
mg_fit_quantileSkewedStable <- function( x )
{
	#return( mg_fit_quantileStable( x ) );
	return( mg_dists_stable.fit.quantile( x ) );
}

## Fits Symmetric Stable distributions with McCulloch's Quantile method.
mg_fit_quantileSymStable <- function( x )
{
	#res <- mg_fit_quantileStable( x );
	res <- mg_dists_stable.fit.quantile( x );

	# removes unused information
	res$beta <- res$gamma <- res$delta <- NULL;

	return( res );
}

## Fits Stable distributions with McCulloch's Quantile method.
#mg_fit_quantileStable <- function( x )
#{
#	return( mg_dists_stable.fit.quantile( x, doplot = FALSE ) );
#}

## Fits Student's t distributions with MLE method.
mg_fit_mleStudentT <- function( x )
{
	library("MASS");

	if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
	{
		stop("[MG::FIT::MLESTUDENTT] 'x' must be a non-empty numeric vector");
	}

	x <- MASS::fitdistr( x, "t" );

	return(
		list(
			df = round( as.numeric( x$estimate["df"] ) )
		)
	);
}

mg_fit_mleUniform <- function( x )
{
	library( MASS );

	#dots <- names(list(...))
	#dots <- dots[!is.element(dots, c("upper", "lower"))]
	if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
	{
		stop("[MG::FIT::MLEUNIFORM] 'x' must be a non-empty numeric vector");
	}

	x <- MASS::fitdistr( x, "uniform" ); #FIXME: not implemented

	return(
		list(
			min = as.numeric( x$estimate["min"] ),
			max = as.numeric( x$estimate["max"] )
		)
	);
}

mg_fit_mleWeibull <- function( x )
{
	library( MASS );

	if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
	{
		stop("[MG::FIT::MLEWEIBULL] 'x' must be a non-empty numeric vector");
	}

	# Ensures x is in the range (0,+Inf)
	x <- x[ x > 0 & is.finite(x) ];

	x <- MASS::fitdistr( x, "weibull" );

	return(
		list(
			shape = as.numeric( x$estimate["shape"] ),
			scale = as.numeric( x$estimate["scale"] )
		)
	);
}

## Estimates Weibull parameters through the least-squared method (see [1])
##
## [1] Ross S., "Introduction to Probability and Statistics for Engineering and
##     Scientists, 2nd Ed."
##
mg_fit_lsqWeibull <- function( x )
{
	# Ensures x is in the range (0,+Inf)
	x <- x[ x > 0 & is.finite(x) ];
	x <- x[ order(x) ]; # we use order-statistics
	n <- length(x);
	y <- c(); # y[i] ~ log(invscale)+shape*log(x[i])
	for ( i in 1:n )
	{
		#y[i] <- log( log( (n + 1)/(n-i+1) ) ); # Method 1
		y[i] <- log( sum( 1 / c(n:(n-i+1) ) ) ); # Method 2
	}
	mlogx <- mean( log( x ) );
	my <- mean( y );
	shape <- ( sum( y*log(x) ) - n*my*mlogx ) / ( sum(log(x)^2) - n*(mlogx^2) );
	invscale <- exp( my - shape * mlogx );

	return( list( shape = shape, scale = 1/invscale ) );
}

## MG_FIT_HILLESTIMATOR(x,k=1,conf.level=NULL)
##
## SUMMARY
##  Calculates the Hill estimator.
##
## DESCRIPTIO
##  Let x_1,...,x_N a sample, N=length(x) and k<N, the Hill estimator is:
##   \frac{1}{k}\sum_{i=1}^k \log(\frac{x_{(i)}{x_{(k+1)}})
##  where x_{(1)},...,x_{(N)} is the inverse order statistics, that is:
##   x_{(1)} >= x_{(2]} >= ... >= x_{(N)}
##
## PARAMETERS
##  x: data values
##  k.max: max number of estimators
##  conf: confidence level for a conf% confidence interval.
mg_fit_hillEstimator <- function( x, k.max = 1, conf.level = NA )
{
	#x <- as.numeric(x);
	x <- x[ is.finite(x) & x > 0 ];
	x <- sort(x, decreasing = TRUE );

	k.max <- min( k.max, length(x)-1 );
	if ( k.max < 1 )
	{
		stop( "[MG::FIT::HILLESTIMATOR] Number of iterations (k.max) must be greater than zero" );
	}

	#k <- 1:(k.max+1);
	k1 <- 1:k.max;
	k2 <- 2:(k.max+1);
	x <- c( log( x[k1] ), log( x[k.max+1] ) );
	#inv.alpha.hat <- ( cumsum( x[ -(k.max+1) ] ) - k[ -(k.max+1) ] * x[ k[-1] ] )  / k[ -(k.max+1) ];
	#inv.alpha.hat <- ( cumsum( x[k1] ) - k1 * x[k2] )  / k1;
	inv.alpha.hat <- cumsum( x[k1] ) / k1 - x[k2];
	u <- l <- NA;
	if ( !is.na( conf.level ) )
	{
		# Wald confidence interval (see Haeusler "Assessing Confidence Intervals for the Tail Index by Edgeworth Expansions for the Hill Estimator" (2005))

		a <- 1 - conf.level;
		z <- qnorm( a/2, lower.tail = FALSE );
		#z <- qt( 1 - a/2, df = k[-1]  );
		s <- z * inv.alpha.hat / sqrt(k1);
		u <- inv.alpha.hat + s;
		l <- inv.alpha.hat - s;
		#u <- 1/u;
		#l <- 1/l;
	}
	#alpha.hat <- 1 / inv.alpha.hat;

	#alpha.hat <- (1/kmax) * sum( ( x[k] ) ) - kmax*x[kmax+1];
	#print( alpha.hat );
	return( list( statistic = inv.alpha.hat, k.max = k.max, conf.int = list(level = conf.level, upper = u, lower = l ) ) );
}

#mg_fit_hillEstimatorOld <- function( x, k = NULL )
#{
#	if ( missing( x ) || length( x ) == 0 || mode( x ) != "numeric" ) 
#	{
#		stop( "'x' must be a non-empty numeric vector" );
#	}
#
#	x <- x[ x > 0 & is.finite(x) ]; # checks for zero value at denominator
#	x <- sort( x ); # we work with order-statistics (x[1] <= x[2] <= ... )
#
#	n <- length( x );
#
#	if ( missing( k ) || is.null( k ) )
#	{
#		k <- n-1;
#	}
#	if ( k <= 0 || k >= n )
#	{
#		stop( "'k' must be  greater than 0 and less than sample ('x') size." );
#	}
#
##	# checks for zero value at denominator
##	#if ( x[n-k+1] == 0 )
##	kk <- k;
##	while ( kk >= 1 )
##	{
##		nn <- n;
##		while ( nn > kk && x[nn-kk] == 0 )
##		{
##			# ignore
##			nn <- nn - 1;
##		}
##		if ( nn > kk )
##		{
##			# found the right k
##
##			if ( nn != n )
##			{
##				warning( c( "Parameter 'n' has been resized to ", nn ) );
##			}
##			n <- nn;
##			break;
##		}
##		kk <- kk - 1;
##	}
##	if ( kk == 0 )
##	{
##		stop( "Too many zeros in data vector." );
##	}
##	if ( kk != k )
##	{
##		warning( c( "Parameter 'k' has been resized to ", kk ) );
##	}
#
#	s <- 0;
#	for ( i in 1:k )
#	{
#		# checks for zero value at numerator
#		if ( x[n-i+1] == 0 )
#		{
#			# ignore
#			next;
#		}
#
#		# checks for negative values in log operand
#		#if ( x[n-i+1] < 0 || x[n-k+1] < 0 )
#		if ( x[n-i+1] < 0 || x[n-k] < 0 )
#		{
#			# ignore
#			next;
#		}
#		#s <- s + log( x[n-i+1]/x[n-h+1] );
#		s <- s + log( x[n-i+1]/x[n-k] );
#	}
#
#	# s/k is the extreme value (gamma)index estimator
#	# To obtain the tail index, simply return k/s
#
#	#return( c( shape = k / s ) );
#	return( list( shape = k / s, xmin = min( x ) ) ); #FIXME: variante min( x[ n - h + 1 ] )
#}

## Fit a power-law distribution as a Yule process (a pure birth process).
## See "Power laws, Pareto distributions and Zipf's law" [Newman].
## Params
## * x: data
## * xmin: minimum value
## * start: starting value for the exponent (shape parameter)
## * ...: additional parameters to MLE function.
mg_fit_yuleEstimator <- function( x, xmin = NULL, start = 2, ... )
{
	if (length(x) == 0)
	{
		stop("[MG::FIT::YULEESTIMATOR] zero length vector")
	}
	if (length(x) == 1) {
		stop("[MG::FIT::YULEESTIMATOR] vector should be at least of length two")
	}

	require( "stats4" );

	x <- x[ x > 0 & is.finite(x) ]; # remove zero values

	if (is.null(xmin)) { xmin <- min(x) }

	n <- length(x)
	x <- x[ x >= xmin ]
	if (length(x) != n) {
		warning("[MG::FIT::YULEESTIMATOR] too small values eliminated from vector")
		n <- length(x)
	}

	#  mlogl <- function(alpha) {
	#    if (xmin > 1) {
	#      C <- 1/(1/(alpha-1)-sum(beta(1:(xmin-1), alpha)))
	#    } else {
	#      C <- alpha-1
	#    }
	#    -n*log(C)-sum(lbeta(x, alpha))
	#  }

	mlogl <- function( alpha ) {
		C <- 1 / sum( ( xmin:10000 )^( -alpha ) );
		return(- n * log( C ) + alpha * sum( log( x ) ) );
	};

	mleres <- stats4::mle( mlogl, start = list( alpha = start ), ... );

	return(
		list(
			shape = as.numeric( (mleres@coef)[["alpha"]] ),
			xmin = xmin
		)
	);
}

#mg_fit_multifit <- function( x, distrs, methods, plotType = NULL, title = NULL, xlabel = NULL, ylabel = NULL, colors = NULL )
mg_fit_multifit <- function( x, distrsMethods )
{
	# preconditions
	if ( missing( distrsMethods ) )
	{
		stop( "[MG::FIT::MULTIFIT] You must specify at least one distribution and one method for each distribution (parameter 'distrsMethods')." );
	}
	distrs <- names( distrsMethods );

	x <- na.omit(x); # remove NA (can cause problem on fitting methods) FIXME
	estimates <- NULL
	for ( d in distrs )
	{
		estimates[[ d ]] <- NULL;
		#estimates[[ d ]] <- list();

		for ( m in distrsMethods[[ d ]] )
		{
			est <- try(
				mg_fit_method( d, m, x )
			);
#print("AIOOO");
#print(d);
#print(m);
#print(est);
			#if ( !is( est, "try-error" ) && !any( is.na( est ) ) )
			if ( !is( est, "try-error" ) && all( is.finite( unlist(est) ) ) )
			{
				estimates[[ d ]] <- list( methods = c( m ), estimate = est );
				break;
			} else {
				msg <- "";
				if ( is( est, "try-error" ) )
				{
					msg <- geterrmessage();
				} else {
					params.str <- paste( est, collapse="," );
					msg <- paste( "Non-Finite Number(s) returned (", params.str, ")", sep="");
				}
				warning( "[MG::FIT::MULTIFIT] Failed to fit parameters for distribution '", d, "' and method '", m, "': ", msg );
			}
		}
		if ( is.null( estimates[[ d ]] ) || length( estimates[[ d ]] ) == 0 )
		{
			warning( "[MG::FIT::MULTIFIT] Failed to fit parameters for distribution '", d, "'." );
		}
	}

	return( estimates );
}
