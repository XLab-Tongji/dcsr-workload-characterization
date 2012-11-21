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

## MG_DISTS
##
## SUMMARY
##  Additional probability distributions.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

source( "lib/mg/mg_consts.R" );
source( "lib/mg/mg_dists_cph.R" ); # Continous Phase-Type
source( "lib/mg/mg_dists_frechet.R" ); # Frechet (Type II Extreme Value)
source( "lib/mg/mg_dists_gev.R" ); # Generalized Extreme Value
source( "lib/mg/mg_dists_gpd.R" ); # Generalized Pareto
source( "lib/mg/mg_dists_gumbel.R" ); # Gumbel (Type I Extreme Value)
source( "lib/mg/mg_dists_hyperexp.R" ); # Hyper-Exponential
source( "lib/mg/mg_dists_pareto.R" ); # Pareto
source( "lib/mg/mg_dists_rweibull.R" ); # Reverse Weibull (Type III Extreme Value)
source( "lib/mg/mg_dists_stable.R" ); # Skewed and Symmetric Stable

####
## Anderson-Darling
####

## Anderson-Darling CDF
##
## [1] Marsaglia, G; Marsaglia JCW; (2004) "Evaluating the Anderson Darling
##
#mg_dists_pad <- function( x, n )
#{
#	p <- mg_dists_adinf( x );
#	p <- p + mg_dists_aderrfix( p, n );
#
#	return( p );
#}
#
### The "ADinf" function see [1].
###
### [1] Marsaglia, G; Marsaglia JCW; (2004) "Evaluating the Anderson Darling
###
#mg_dists_adinf <- function( x )
#{
#	y <- array( dim = length(x) );
#	idx1 <- which( x < 2 );
#	idx2 <- which( x >= 2 );
#
#	if ( length( idx1 ) )
#	{
#		x1 <- x[ idx1 ];
#		# max |error| < .000002 for x<2, (p=.90816...)
#		y[ idx1 ] <- exp( -1.2337141 / x1 ) / sqrt( x1 ) * ( 2.00012 + ( .247105 - ( .0649821 - ( .0347962 - ( .011672 - .00168691 * x1 ) * x1 ) * x1 ) * x1 ) * x1 );
#	}
#	if ( length( idx2 ) )
#	{
#		x2 <- x[ idx2 ];
#		# max |error|<.0000008 for 4<x<infinity
# 		y[ idx2 ] <- exp( -exp( 1.0776 - ( 2.30695 - ( .43424 - ( .082433 - ( .008056 - .0003146 * x2 ) * x2 ) * x2 ) * x2 ) * x2 ) );
#	}
#
#	return( y );
#}
#
### The "errfix" function see [1].
###
### [1] Marsaglia, G; Marsaglia JCW; (2004) "Evaluating the Anderson Darling
###
#mg_dists_aderrfix <- function( x, n )
#{
#	if ( is.vector( x ) )
#	{
#		if ( is.numeric( n ) )
#		{
#			n <- seq( n, n, length.out = length(x) );
#		} else if ( length( n ) != length( x ) ) {
#			rlen <- length( x ) - length( n );
#			if ( rlen > 0 )
#			{
#				n <- c( n, seq( n[1], n[1], length.out = rlen ) );
#			} else {
#				n <- n[ 1 : length(x) ];
#			}
#		}
#	}
#
#	idx1 <- which( x >= 0.8 );
#	idx2 <- which( x < 0.8 );
#
#	if ( length( idx1 ) > 0 )
#	{
#		y <- x[ idx1 ];
#		y <- ( ( -130.2137 + ( 745.2337 - ( 1705.091 - ( 1950.646 - ( 1116.360 - 255.7844 * y ) * y ) * y ) * y ) * y ) / n[ idx1 ] );
#		x[ idx1 ] <- y;
#	}
#
#	if ( length( idx2 ) > 0 )
#	{
#		c <- 0.1265 + 0.1757 / n[ idx2 ];
#		y <- x[ idx2 ];
#
#		idx3 <- which( y < c );
#		idx4 <- which( y >= c );
#
#		if ( length( idx3 ) > 0 )
#		{
#			n3 <- n[ idx3 ];
#			t <- y[ idx3 ] / c;
#			t <- sqrt( t ) * ( 1 - t ) * ( 49 * t - 102 );
#			y[ idx3 ] <- t * ( 0.0037 / (n3 * n3) + 0.00078 / n3 + 0.00006 ) / n3;
#		}
#		if ( length( idx4 ) > 0 )
#		{
#			n4 <- n[ idx4 ];
#			t <- ( y[ idx4 ] - c ) / ( 0.8 - c );
#			t <- sqrt( t ) * ( 1 - t ) * ( 49 * t - 102 );
#			t <- -0.00022633 + ( 6.54034 - ( 14.6538 - ( 14.458 - ( 8.259 - 1.91864 * t ) * t ) * t ) * t ) * t;
#			y[ idx4 ] <- t * ( 0.04213 + 0.01365 / n4 ) / n4;
#		}
#		x[ idx2 ] <- y;
#	}
#
#	return( x );
#}

####
## Utility functions.
####

## MG_DISTS_PDFNAME(dist)
##
## SUMMARY
##  Returns the probability distribution function name of a given distribution
##  name.
##
## PARAMETERS
##  * dist: the probability distribution name; see 'mg_consts.R' for possible
##    probability distribution names.
##
## RETURN
##  The name of probability distribution function.
##
## DESCRIPTION
##  Given the name 'dist' of a probability distribution, returns the related
##  probability distribution function name. E.g. for 'normal' distribution
##  (MG_CONSTS_NORMAL_DIST)j returns 'dnorm'.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##
mg_dists_pdfName <- function( dist )
{
	if ( dist == MG_CONSTS_BETA_DIST )
	{
		return( "dbeta" );
	}
	if ( dist == MG_CONSTS_BINOMIAL_DIST )
	{
		return( "dbinom" );
	}
	if ( dist == MG_CONSTS_CAUCHY_DIST )
	{
		return( "dcauchy" );
	}
	if ( dist == MG_CONSTS_CHISQUARED_DIST )
	{
		return( "dchisq" );
	}
	if ( dist == MG_CONSTS_EXPONENTIAL_DIST )
	{
		return( "dexp" );
	}
	if ( dist == MG_CONSTS_FISHERF_DIST )
	{
		return( "df" );
	}
	if ( dist == MG_CONSTS_FRECHET_DIST )
	{
		return( "mg_dists_dfrechet" );
	}
	if ( dist == MG_CONSTS_GAMMA_DIST )
	{
		return( "dgamma" );
	}
	if ( dist == MG_CONSTS_GENERALIZEDEXTREMEVALUE_DIST )
	{
		return( "mg_dists_dgev" );
	}
	if ( dist == MG_CONSTS_GENERALIZEDPARETO_DIST )
	{
		return( "mg_dists_dgpd" );
	}
	if ( dist == MG_CONSTS_GUMBEL_DIST )
	{
		return( "mg_dists_dgumbel" );
	}
	if ( dist == MG_CONSTS_HYPEREXPONENTIAL_DIST )
	{
		return( "mg_dists_dhyperexp" );
	}
	if ( dist == MG_CONSTS_LOGISTIC_DIST )
	{
		return( "dlogis" );
	}
	if ( dist == MG_CONSTS_LOGNORMAL_DIST )
	{
		return( "dlnorm" );
	}
	if ( dist == MG_CONSTS_NORMAL_DIST )
	{
		return( "dnorm" );
	}
	if ( dist == MG_CONSTS_PARETO_DIST )
	{
		return( "mg_dists_dpareto" );
	}
	if ( dist == MG_CONSTS_CONTINUOUSPHASETYPE_DIST )
	{
		return( "mg_dists_dcph" );
	}
	if ( dist == MG_CONSTS_POISSON_DIST )
	{
		return( "dpois" );
	}
	if ( dist == MG_CONSTS_REVERSEWEIBULL_DIST )
	{
		return( "mg_dists_drweibull" );
	}
	if ( dist == MG_CONSTS_SKEWEDSTABLE_DIST )
	{
		return( "mg_dists_dstable" );
	}
	if ( dist == MG_CONSTS_STABLE_DIST )
	{
		return( "mg_dists_dstable" );
	}
	if ( dist == MG_CONSTS_STUDENTT_DIST )
	{
		return( "dt" );
	}
	if ( dist == MG_CONSTS_SYMSTABLE_DIST )
	{
		return( "mg_dists_dsymstb" );
	}
	if ( dist == MG_CONSTS_UNIFORM_DIST )
	{
		return( "duniform" );
	}
	if ( dist == MG_CONSTS_WEIBULL_DIST )
	{
		return( "dweibull" );
	}
	if ( dist == MG_CONSTS_WEIBULL3_DIST )
	{
		return( "mg_dists_dweibull3" );
	}
	return( "" );
}

## MG_DISTS_DDF(dist,x,...)
##
## SUMMARY
##  Returns the value of probability distribution function evaluated at q.
##
## PARAMETERS
##  * dist: the probability distribution name; see @c mg_consts.R for possible
##    probability distribution names.
##  * q: quantile values.
##  * ...: list of distribution parameters.
##
## RETURN
##  The value of probability distribution function evaluated at x.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##
mg_dists_pdf <- function( dist, ... )
{
	pdfname <- mg_dists_pdfName( dist );
	return( do.call( pdfname, list( ... ) ) );
}

## MG_DISTS_CDFNAME(dist)
##
## SUMMARY
##  Returns the cumulative distribution function name of a given distribution
##  name.
##
## PARAMETERS
##  * dist: the cumulative distribution name; see 'mg_consts.R' for possible
##    probability distribution names.
##
## RETURN
##  The name of cumulative distribution function.
##
## DESCRIPTION
##  Given the name 'dist' of a probability distribution, returns the related
##  cumulative distribution function name. E.g. for 'normal' distribution
##  (MG_CONSTS_NORMAL_DIST)j returns 'pnorm'.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##
mg_dists_cdfName <- function( dist )
{
	if ( dist == MG_CONSTS_BETA_DIST )
	{
		return( "pbeta" );
	}
	if ( dist == MG_CONSTS_BINOMIAL_DIST )
	{
		return( "pbinom" );
	}
	if ( dist == MG_CONSTS_CAUCHY_DIST )
	{
		return( "pcauchy" );
	}
	if ( dist == MG_CONSTS_CHISQUARED_DIST )
	{
		return( "pchisq" );
	}
	if ( dist == MG_CONSTS_EXPONENTIAL_DIST )
	{
		return( "pexp" );
	}
	if ( dist == MG_CONSTS_FISHERF_DIST )
	{
		return( "pf" );
	}
	if ( dist == MG_CONSTS_FRECHET_DIST )
	{
		return( "mg_dists_pfrechet" );
	}
	if ( dist == MG_CONSTS_GAMMA_DIST )
	{
		return( "pgamma" );
	}
	if ( dist == MG_CONSTS_GENERALIZEDEXTREMEVALUE_DIST )
	{
		return( "mg_dists_pgev" );
	}
	if ( dist == MG_CONSTS_GENERALIZEDPARETO_DIST )
	{
		return( "mg_dists_pgpd" );
	}
	if ( dist == MG_CONSTS_GUMBEL_DIST )
	{
		return( "mg_dists_pgumbel" );
	}
	if ( dist == MG_CONSTS_HYPEREXPONENTIAL_DIST )
	{
		return( "mg_dists_phyperexp" );
	}
	if ( dist == MG_CONSTS_LOGISTIC_DIST )
	{
		return( "plogis" );
	}
	if ( dist == MG_CONSTS_LOGNORMAL_DIST )
	{
		return( "plnorm" );
	}
	if ( dist == MG_CONSTS_NORMAL_DIST )
	{
		return( "pnorm" );
	}
	if ( dist == MG_CONSTS_PARETO_DIST )
	{
		return( "mg_dists_ppareto" );
	}
	if ( dist == MG_CONSTS_CONTINUOUSPHASETYPE_DIST )
	{
		return( "mg_dists_pcph" );
	}
	if ( dist == MG_CONSTS_POISSON_DIST )
	{
		return( "ppois" );
	}
	if ( dist == MG_CONSTS_REVERSEWEIBULL_DIST )
	{
		return( "mg_dists_prweibull" );
	}
	if ( dist == MG_CONSTS_SKEWEDSTABLE_DIST )
	{
		return( "mg_dists_pstable" );
	}
	if ( dist == MG_CONSTS_STABLE_DIST )
	{
		return( "mg_dists_pstable" );
	}
	if ( dist == MG_CONSTS_STUDENTT_DIST )
	{
		return( "pt" );
	}
	if ( dist == MG_CONSTS_SYMSTABLE_DIST )
	{
		return( "mg_dists_psymstb" );
	}
	if ( dist == MG_CONSTS_UNIFORM_DIST )
	{
		return( "puniform" );
	}
	if ( dist == MG_CONSTS_WEIBULL_DIST )
	{
		return( "pweibull" );
	}
	if ( dist == MG_CONSTS_WEIBULL3_DIST )
	{
		return( "mg_dists_pweibull3" );
	}
	return( "" );
}

## MG_DISTS_CDF(dist,q,...)
##
## SUMMARY
##  Returns the value of cumulative distribution function evaluated at q.
##
## PARAMETERS
##  * dist: the probability distribution name; see @c mg_consts.R for possible
##    probability distribution names.
##  * q: quantile values.
##  * ...: list of distribution parameters.
##
## RETURN
##  The value of cumulative distribution function evaluated at q.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##
mg_dists_cdf <- function( dist, ... )
{
	cdfname <- mg_dists_cdfName( dist );
	return( do.call( cdfname, list( ... ) ) );
}

## MG_DISTS_RVGNAME(dist)
##
## SUMMARY
##  Returns the name of random variate generator function for a given
##  distribution name.
##
## PARAMETERS
##  * dist: the probability distribution name; see 'mg_consts.R' for possible
##    probability distribution names.
##
## RETURN
##  The name of random variate generator function.
##
## DESCRIPTION
##  Given the name 'dist' of a probability distribution, returns the related
##  random variate generator function name. E.g. for 'normal' distribution
##  (MG_CONSTS_NORMAL_DIST)j returns 'rnorm'.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##
mg_dists_rvgName <- function( dist )
{
	if ( dist == MG_CONSTS_BETA_DIST )
	{
		return( "rbeta" );
	}
	if ( dist == MG_CONSTS_BINOMIAL_DIST )
	{
		return( "rbinom" );
	}
	if ( dist == MG_CONSTS_CAUCHY_DIST )
	{
		return( "rcauchy" );
	}
	if ( dist == MG_CONSTS_CHISQUARED_DIST )
	{
		return( "rchisq" );
	}
	if ( dist == MG_CONSTS_EXPONENTIAL_DIST )
	{
		return( "rexp" );
	}
	if ( dist == MG_CONSTS_FISHERF_DIST )
	{
		return( "rf" );
	}
	if ( dist == MG_CONSTS_FRECHET_DIST )
	{
		return( "mg_dists_rfrechet" );
	}
	if ( dist == MG_CONSTS_GAMMA_DIST )
	{
		return( "rgamma" );
	}
	if ( dist == MG_CONSTS_GENERALIZEDEXTREMEVALUE_DIST )
	{
		return( "mg_dists_rgev" );
	}
	if ( dist == MG_CONSTS_GENERALIZEDPARETO_DIST )
	{
		return( "mg_dists_rgpd" );
	}
	if ( dist == MG_CONSTS_GUMBEL_DIST )
	{
		return( "mg_dists_rgumbel" );
	}
	if ( dist == MG_CONSTS_HYPEREXPONENTIAL_DIST )
	{
		return( "mg_dists_rhyperexp" );
	}
	if ( dist == MG_CONSTS_LOGISTIC_DIST )
	{
		return( "rlogis" );
	}
	if ( dist == MG_CONSTS_LOGNORMAL_DIST )
	{
		return( "rlnorm" );
	}
	if ( dist == MG_CONSTS_NORMAL_DIST )
	{
		return( "rnorm" );
	}
	if ( dist == MG_CONSTS_PARETO_DIST )
	{
		return( "mg_dists_rpareto" );
	}
	if ( dist == MG_CONSTS_CONTINUOUSPHASETYPE_DIST )
	{
		return( "mg_dists_rcph" );
	}
	if ( dist == MG_CONSTS_POISSON_DIST )
	{
		return( "rpois" );
	}
	if ( dist == MG_CONSTS_REVERSEWEIBULL_DIST )
	{
		return( "mg_dists_rrweibull" );
	}
	if ( dist == MG_CONSTS_SKEWEDSTABLE_DIST )
	{
		return( "mg_dists_rstable" );
	}
	if ( dist == MG_CONSTS_STABLE_DIST )
	{
		return( "mg_dists_rstable" );
	}
	if ( dist == MG_CONSTS_STUDENTT_DIST )
	{
		return( "rt" );
	}
	if ( dist == MG_CONSTS_SYMSTABLE_DIST )
	{
		return( "mg_dists_rsymstb" );
	}
	if ( dist == MG_CONSTS_UNIFORM_DIST )
	{
		return( "runiform" );
	}
	if ( dist == MG_CONSTS_WEIBULL_DIST )
	{
		return( "rweibull" );
	}
	if ( dist == MG_CONSTS_WEIBULL3_DIST )
	{
		return( "mg_dists_rweibull3" );
	}
	return( "" );
}

## MG_DISTS_RVG(dist,n,...)
##
## SUMMARY
##  Returns the random variate values for the given distribution.
##
## PARAMETERS
##  * dist: the probability distribution name; see 'mg_consts.R' for possible
##    probability distribution names.
##  * n: number of variates to be generated.
##  * ...: list of distribution parameters.
##
## RETURN
##  The values of random variates.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##
mg_dists_rvg <- function( dist, ... )
{
	rvgname <- mg_dists_rvgName( dist );
	return( do.call( rvgname, list( ... ) ) );
}

## MG_DISTS_INVCDFNAME(dist)
##
## SUMMARY
##  Returns the name of inverse CDF (a.k.a. quantile function) for a given
##  distribution name.
##
## PARAMETERS
##  * dist: the probability distribution name; see 'mg_consts.R' for possible
##    probability distribution names.
##
## RETURN
##  The name of inverse CDF.
##
## DESCRIPTION
##  Given the name 'dist' of a probability distribution, returns the related
##  inverse CDF name. E.g. for 'normal' distribution (MG_CONSTS_NORMAL_DIST)j
##  returns 'qnorm'.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##
mg_dists_invcdfName <- function( dist )
{
	if ( dist == MG_CONSTS_BETA_DIST )
	{
		return( "qbeta" );
	}
	if ( dist == MG_CONSTS_BINOMIAL_DIST )
	{
		return( "qbinom" );
	}
	if ( dist == MG_CONSTS_CAUCHY_DIST )
	{
		return( "qcauchy" );
	}
	if ( dist == MG_CONSTS_CHISQUARED_DIST )
	{
		return( "qchisq" );
	}
	if ( dist == MG_CONSTS_EXPONENTIAL_DIST )
	{
		return( "qexp" );
	}
	if ( dist == MG_CONSTS_FISHERF_DIST )
	{
		return( "qf" );
	}
	if ( dist == MG_CONSTS_FRECHET_DIST )
	{
		return( "mg_dists_qfrechet" );
	}
	if ( dist == MG_CONSTS_GAMMA_DIST )
	{
		return( "qgamma" );
	}
	if ( dist == MG_CONSTS_GENERALIZEDEXTREMEVALUE_DIST )
	{
		return( "mg_dists_qgev" );
	}
	if ( dist == MG_CONSTS_GENERALIZEDPARETO_DIST )
	{
		return( "mg_dists_qgpd" );
	}
	if ( dist == MG_CONSTS_GUMBEL_DIST )
	{
		return( "mg_dists_qgumbel" );
	}
	if ( dist == MG_CONSTS_HYPEREXPONENTIAL_DIST )
	{
		return( "mg_dists_qhyperexp" );
	}
	if ( dist == MG_CONSTS_LOGISTIC_DIST )
	{
		return( "qlogis" );
	}
	if ( dist == MG_CONSTS_LOGNORMAL_DIST )
	{
		return( "qlnorm" );
	}
	if ( dist == MG_CONSTS_NORMAL_DIST )
	{
		return( "qnorm" );
	}
	if ( dist == MG_CONSTS_PARETO_DIST )
	{
		return( "mg_dists_qpareto" );
	}
	if ( dist == MG_CONSTS_CONTINUOUSPHASETYPE_DIST )
	{
		return( "mg_dists_qcph" );
	}
	if ( dist == MG_CONSTS_POISSON_DIST )
	{
		return( "qpois" );
	}
	if ( dist == MG_CONSTS_REVERSEWEIBULL_DIST )
	{
		return( "mg_dists_qrweibull" );
	}
	if ( dist == MG_CONSTS_SKEWEDSTABLE_DIST )
	{
		return( "mg_dists_qstable" );
	}
	if ( dist == MG_CONSTS_STABLE_DIST )
	{
		return( "mg_dists_qstable" );
	}
	if ( dist == MG_CONSTS_STUDENTT_DIST )
	{
		return( "qt" );
	}
	if ( dist == MG_CONSTS_SYMSTABLE_DIST )
	{
		return( "mg_dists_qsymstb" );
	}
	if ( dist == MG_CONSTS_UNIFORM_DIST )
	{
		return( "quniform" );
	}
	if ( dist == MG_CONSTS_WEIBULL_DIST )
	{
		return( "qweibull" );
	}
	if ( dist == MG_CONSTS_WEIBULL3_DIST )
	{
		return( "mg_dists_qweibull3" );
	}
	return( "" );
}

mg_dists_invcdf <- function( dist, ... )
{
	invcdfname <- mg_dists_invcdfName( dist );
	return( do.call( invcdfname, list( ... ) ) );
}

mg_dists_toString <- function( dist )
{
	if ( dist == MG_CONSTS_BETA_DIST )
	{
		return( "Beta" );
	}
	if ( dist == MG_CONSTS_BINOMIAL_DIST )
	{
		return( "Binomial" );
	}
	if ( dist == MG_CONSTS_CAUCHY_DIST )
	{
		return( "Cauchy" );
	}
	if ( dist == MG_CONSTS_CHISQUARED_DIST )
	{
		return( "Chi-Square" );
	}
	if ( dist == MG_CONSTS_EXPONENTIAL_DIST )
	{
		return( "Exp" );
	}
	if ( dist == MG_CONSTS_FISHERF_DIST )
	{
		return( "Fisher's F" );
	}
	if ( dist == MG_CONSTS_FRECHET_DIST )
	{
		return( "Frechet" );
	}
	if ( dist == MG_CONSTS_GAMMA_DIST )
	{
		return( "Gamma" );
	}
	if ( dist == MG_CONSTS_GENERALIZEDEXTREMEVALUE_DIST )
	{
		return( "GEV" );
	}
	if ( dist == MG_CONSTS_GENERALIZEDPARETO_DIST )
	{
		return( "GPD" );
	}
	if ( dist == MG_CONSTS_GUMBEL_DIST )
	{
		return( "Gumbel" );
	}
	if ( dist == MG_CONSTS_HYPEREXPONENTIAL_DIST )
	{
		return( "Hyper-Exp" );
	}
	if ( dist == MG_CONSTS_LOGISTIC_DIST )
	{
		return( "Logistic" );
	}
	if ( dist == MG_CONSTS_LOGNORMAL_DIST )
	{
		return( "Log-Normal" );
	}
	if ( dist == MG_CONSTS_NORMAL_DIST )
	{
		return( "Normal" );
	}
	if ( dist == MG_CONSTS_PARETO_DIST )
	{
		return( "Pareto" );
	}
	if ( dist == MG_CONSTS_CONTINUOUSPHASETYPE_DIST )
	{
		return( "CPH" );
	}
	if ( dist == MG_CONSTS_POISSON_DIST )
	{
		return( "Poisson" );
	}
	if ( dist == MG_CONSTS_REVERSEWEIBULL_DIST )
	{
		return( "Reverse-Weibull" );
	}
	if ( dist == MG_CONSTS_SKEWEDSTABLE_DIST )
	{
		return( "Skewed-Stable" );
	}
	if ( dist == MG_CONSTS_STABLE_DIST )
	{
		return( "Stable" );
	}
	if ( dist == MG_CONSTS_STUDENTT_DIST )
	{
		return( "Student's t" );
	}
	if ( dist == MG_CONSTS_SYMSTABLE_DIST )
	{
		return( "Symmetric-Stable" );
	}
	if ( dist == MG_CONSTS_UNIFORM_DIST )
	{
		return( "Uniform" );
	}
	if ( dist == MG_CONSTS_WEIBULL_DIST )
	{
		return( "Weibull" );
	}
	if ( dist == MG_CONSTS_WEIBULL3_DIST )
	{
		return( "Weibull-3" );
	}
	return( "" );
}

mg_dists_isDiscrete <- function( dist )
{
	if (
		dist == MG_CONSTS_BINOMIAL_DIST
		|| dist == MG_CONSTS_POISSON_DIST
#		|| dist == MG_CONSTS_DPH_DIST
	)
	{
		return( TRUE );
	}
	return( FALSE );
}

############## EXPERIMENTAL SECTION ##########################################

####
## Truncated Weibull
####

mg_dists_dweibull.trunc <- function( x, shape, scale=1, trunc. = Inf, log = FALSE )
{
	ln.dens <- dweibull( x, shape, scale, log = TRUE ) - pweibull( trunc., shape, scale = 1, lower.tail = TRUE, log.p = TRUE )

	if ( any( oops <- ( x > trunc. ) ) )
	{
		ln.dens[ oops ] <- -Inf
	}
	if ( log )
	{
		return( ln.dens )
	}

	return( exp( ln.dens ) )

#FIXME
#If you want to estimate the truncation point, that will be a more
#difficult problem. For that, I suggest you compute the max of your data
#and parameterize the truncated density with a parameter like
#"log.trunc.over.max" so "trunc." in the above example is computed as
#(max+exp(log.trunc.over.max)).
}


