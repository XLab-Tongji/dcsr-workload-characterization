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

## MG_GOF_KS
##
## SUMMARY
##  Kolmogorov-Smirnov Goodness-of-Fit tests.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

source( "lib/mg/mg_dists.R" );
source( "lib/mg/mg_fit.R" );

## MG_GOF_KS.TEST(x,dist,...,alternative=c("two.sided","less","greater"),exact=NULL)
##
## SUMMARY
##  1-Sample Kolmogorov-Smirnov GoF test.
##
## PARAMS
##  * x: data values.
##  * dist: theoretical distribution name.
##  * ...: parameters to pass to the diven distribution @c dist.
##  * alternative: type of alternative hypothesis.
##  * exact: NULL or logical value indicating whether to compute an exact
##    p-value.
##
## RETURN
##  A list with the following fields:
##  * statistic: the KS test statistic.
##  * p.value: the p-value of the test.
##  * alternative: a string describing the type of alternative hypothesis.
##  * method: a string indicating what type of test was performed.
##  * data.name: a string giving the name(s) of the data.
##
## SEE ALSO
##  ks.test, mg_gof_ks.test.boot, mg_gof_ks.test.sim, mg_gof_ks.test2, mg_gof_ad
##
mg_gof_ks.test <- function( x, dist, ..., alternative = c("two.sided", "less", "greater"), exact = NULL )
{
	#return( ks.test( x, mg_dists_cdfName( dist ), ..., alt = alternative, ex = exact ) );
	return( ks.test( x, mg_dists_cdfName( dist ), ..., alternative = alternative, exact = exact ) );
}

## MG_GOF_KS.TEST.BOOT.PAR(x,dist,...,n.boot=999,alternative=c("two.sided","less","greater"),exact=NULL)
##
## SUMMARY
##  1-Sample Kolmogorov-Smirnov GoF test with Parametric Bootstraping.
##
## PARAMS
##  * x: data values.
##  * dist: theoretical distribution name.
##  * ...: parameters to pass to the diven distribution @c dist.
##  * n.boot: number of bootstrap iterations.
##  * alternative: type of alternative hypothesis.
##  * exact: NULL or logical value indicating whether to compute an exact
##    p-value.
##
## RETURN
##  A list with the following fields:
##  * statistic: the KS test statistic without bootstrapping.
##  * p.value: the p-value of the test without bootstrapping.
##  * p.value.boot: the p-value of the test calcolated with bootstrap method.
##  * alternative: a string describing the type of alternative hypothesis.
##  * method: a string indicating what type of test was performed.
##  * data.name: a string giving the name(s) of the data.
##
## SEE ALSO
##  ks.test, mg_gof_ks.test, mg_gof_ks.test.sim, mg_gof_ks.test2, mg_gof_ad
##
mg_gof_ks.test.boot.par <- function( x, dist, ..., n.boot = 999, alternative = c("two.sided", "less", "greater"), exact = NULL, fit.methods = NULL )
{
	n.x <- length(x);

	# Estimates parameters for dist from X^* with the same method used
	# for parameters passed as argument
	if ( is.null( fit.methods ) || length( fit.methods ) == 0 )
	{
		fit.methods <- mg_fit_methodsNameForDistrs( dist );
	} else {
		tmp <- fit.methods;
		fit.methods <- list();
		fit.methods[[ dist ]] <- tmp;
	}

	# Perform a K-S test from X and distr (whose parameters are estimated from X)
	cdf <- mg_dists_cdfName( dist );
	ks <- try( ks.test( x, cdf, ..., alternative = alternative, exact = exact ) );
	if ( is( ks, "try-error" ) )
	{
		stop( "[MG::GOF::KS::TEST::BOOT::PAR] Error while performing K-S test for distribution '", dist, "' on original sample: ", geterrmessage() );
	}


	ks.pval <- 0;
	#ks.list <- c();
	boot.mean <- 0;
	boot.var <- 0;
	for ( i in 1:n.boot )
	{
		# Samples from the theoretical distribution (=> the bootstrap sample X^*)
		x.star <- try( mg_dists_rvg( dist, n.x, ... ) );
		if ( is( x.star, "try-error" ) )
		{
			warning( "[MG::GOF::KS.TEST.BOOT.PAR] Error while sampling from distribution '", dist, "' in bootstrap iteration #", i, ": ", geterrmessage(), "." );
			next;
		}

#		# Estimates parameters for dist from X^* with the same method used
#		# for parameters passed as argument
#		if ( is.null( fit.methods ) || length( fit.methods ) == 0 )
#		{
#			fit.methods <- mg_fit_methodsNameForDistrs( dist );
#		}
#		else
#		{
#			tmp <- fit.methods;
#			fit.methods <- list();
#			fit.methods[[ dist ]] <- tmp;
#		}
		theta.star <- try( mg_fit_multifit( x.star, fit.methods ) );
		if ( is( theta.star, "try-error" ) )
		{
			warning( "[MG::GOF::KS.TEST.BOOT.PAR] Error while fitting distribution '", dist, "' in a bootstrap iteration #", i, ": ", geterrmessage(), "." );
			next;
                } else if ( is.null( theta.star ) || length( theta.star ) == 0 ) {
                        warning( "[MG::GOF::KS.TEST.BOOT.PAR] Error while fitting distribution '", dist, "' in bootstrap iteration #", i, ": returned NULL." );
                        next;
                }

		# Performs a K-S test between the bootstrap sample X^* and the distribution dist with parameters theta^*
		#ks.boot <- ks.test( x, x.boot, alternative = alternative, exact = exact );
		ks.star <- try( do.call( "ks.test", c( list(x = x.star, y = cdf, theta.star[[dist]]$estimate , list(alternative = alternative, exact = exact) ) ) ) );
                if ( is( ks.star, "try-error" ) )
                {
                        stop( "[MG::GOF::KS::TEST::BOOT::PAR] Error while performing K-S test for distribution '", dist, "' on ", i, "-th bootstrap sample: ", geterrmessage() );
                }

		#ks.list[i] <- ks.boot$statistic;
		boot.mean <- ks.star$statistic;
		boot.var <- ks.star$statistic^2;

		if ( ks.star$statistic >= ks$statistic )
		{
			ks.pval <- ks.pval + 1;
		}
	}
	#ks.pval <- (1+ks.pval) / (n.boot + 1);
	ks.pval <- ks.pval / (n.boot + 1);
	boot.mean <- boot.mean/n.boot;
	boot.var <- 1/(n.boot-1) * (boot.var - n.boot*boot.mean^2);

	#hist( ks.list)
	#abline( v=c(ks$statistic) );
	return(
		list(
			statistic = ks$statistic,
			p.value = ks$p.value,
			p.value.boot = ks.pval,
			statistic.boot.mean = boot.mean,
			statistic.boot.sd = sqrt(boot.var),
			method = ks$method,
			alternative = ks$alternative,
			data.name = ks$data.name
		)
	);
}

## MG_GOF_KS.TEST.BOOT.NONPAR(x,dist,...,n.boot=999,alternative=c("two.sided","less","greater"),exact=NULL)
##
## SUMMARY
##  1-Sample Kolmogorov-Smirnov GoF test with Non-Parametric Bootstrapping.
##
## PARAMS
##  * x: data values.
##  * dist: theoretical distribution name.
##  * ...: parameters to pass to the diven distribution @c dist.
##  * n.boot: number of bootstrap iterations.
##  * alternative: type of alternative hypothesis.
##  * exact: NULL or logical value indicating whether to compute an exact
##    p-value.
##
## RETURN
##  A list with the following fields:
##  * statistic: the KS test statistic without bootstrapping.
##  * p.value: the p-value of the test without bootstrapping.
##  * p.value.boot: the p-value of the test calcolated with bootstrap method.
##  * alternative: a string describing the type of alternative hypothesis.
##  * method: a string indicating what type of test was performed.
##  * data.name: a string giving the name(s) of the data.
##
## SEE ALSO
##  ks.test, mg_gof_ks.test, mg_gof_ks.test.sim, mg_gof_ks.test2, mg_gof_ad
##
mg_gof_ks.test.boot.nonpar <- function( x, dist, ..., n.boot = 999, alternative = c("two.sided", "less", "greater"), exact = NULL )
{
	n.x <- length(x);

	cdf <- mg_dists_cdfName( dist );
	ks <- ks.test( x, cdf, ..., alternative = alternative, exact = exact );

	ks.pval <- 0;
	ks.list <- c();
	for ( i in 1:n.boot )
	{
		# Samples from the theoretical distribution
		y <- mg_dists_rvg( dist, n.x, ... );

		# Performs a K-S test
		ks.boot <- ks.test( x, y, alternative = alternative, exact = exact );

		ks.list[i] <- ks.boot$statistic;

		if ( ks.boot$statistic >= ks$statistic )
		{
			ks.pval <- ks.pval + 1;
		}
	}
	#ks.pval <- (1+ks.pval) / (n.boot + 1);
	ks.pval <- ks.pval / (n.boot + 1);

	#hist( ks.list)
	#abline( v=c(ks$statistic) );
	return(
		list(
			statistic = ks$statistic,
			p.value = ks$p.value,
			p.value.boot = ks.pval,
			method = ks$method,
			alternative = ks$alternative,
			data.name = ks$data.name
		)
	);
}

## MG_GOF_KS.TEST.SIM(x,dist,...,n.sim=1)
##
## SUMMARY
##  1-Sample simulated Kolmogorov-Smirnov GoF test.
##
## PARAMS
##  * x: the data values
##  * dist: the distribution name agaist which you would test the GoF.
##  * n.sim: number of simulation of 2-Sample KS test.
##
## RETURN
##  A list with the following fields:
##  * statistic: the KS test statistic.
##  * p.value: the p-value of the test.
##  * alternative: a string describing the type of alternative hypothesis.
##  * method: a string indicating what type of test was performed.
##  * data.name: a string giving the name(s) of the data.
##
## DESCRIPTION
##  Performs a 1-Sample Kolmogorov-Smirnov test performing @c n.sim 2-Samples
##  Kolmogorov-Smirnov test, where one sample is the given parameter @c x,
##  and the other is a random variate sample extracted from given distribution
##  @c dist.
##
## SEE ALSO
##  ks.test, mg_gof_ks.test.sim, mg_gof_ks.test2, mg_gof_ad
##
mg_gof_ks.test.sim <- function( x, dist, ..., n.sim = 1000, alternative = c("two.sided", "less", "greater"), exact = NULL )
{
	n <- length(x);
	ks.stat.mu <- 0;
	ks.stat.sd <- 0;
	ks.pval.mu <- 0;
	ks.pval.sd <- 0;
	ks.method <- "";
	ks.alt <- "";
	ks.dname <- "";

	for ( k in 1:n.sim )
	{
		ks <- mg_gof_ks.test2( x, mg_dists_rvg( dist, n, ... ), alt = alternative, ex = exact );

		ks.stat.mu <- ks.stat.mu + ks$statistic;
		ks.stat.sd <- ks.stat.sd + ks$statistic^2;
		ks.pval.mu <- ks.pval.mu + ks$p.value;
		ks.pval.sd <- ks.pval.sd + ks$p.value^2;

		if ( k == 1 )
		{
			ks.method <- ks$method;
			ks.alt <- ks$alternative;
			ks.dname <- ks$data.name;
		}
	}
	ks.stat.mu <- ks.stat.mu/n.sim;
	ks.pval.mu <- ks.pval.mu/n.sim;
	if ( n.sim > 1 )
	{
		ks.stat.sd <- sqrt( ( ks.stat.sd - n.sim*ks.stat.mu^2 ) / (n.sim-1) );
		ks.pval.sd <- sqrt( ( ks.pval.sd - n.sim*ks.pval.mu^2 ) / (n.sim-1) );
	} else {
		ks.stat.sd <- 0;
		ks.pval.sd <- 0;
	}
	names( ks.stat.sd ) <- "";
	names( ks.pval.sd ) <- "";

	#n <- n * n/(n + n);
	#p.value <- 1 - mg_gof_ks.cdf( n, ks.stat.mu );

	return(
		list(
			statistic = ks.stat.mu,
			statistic.sd = ks.stat.sd,
			p.value = ks.pval.mu,
			p.value.sd = ks.pval.sd,
			method = ks.method,
			alternative = ks.alt,
			alternative = alternative,
			data.name = ks.dname
		)
	);
}

## MG_GOF_KS.TEST2(x,y,alternative=c("two-sided","less","greater"),exact=NULL))
##
## SUMMARY
##  2-Samples Kolmogorov-Smirnov GoF test.
##
mg_gof_ks.test2 <- function( x, y, alternative = c("two.sided", "less", "greater"), exact = NULL )
{
	return( ks.test( x, y, alt = alternative, ex = exact ) );
}

## MG_GOF_KS.TEST2.BOOT(x,y,n.boot=100,alternative=c("two-sided","less","greater"),exact=NULL),replace=TRUE)
##
## SUMMARY
##  2-Samples Kolmogorov-Smirnov GoF test with bootstraping.
## 
## PARAMS
##  * x: first data set.
##  * y: second data set.
##  * n.boot: number of resampling.
##  * alternative: type of alternative hypothesis.
##  * exact: @c TRUE if p-value should be calculated exactly.
##  * replace: @c TRUE if resampling should be with replacement.
##
## RETURN
##  A list with the following fields:
##  * statistic: the KS test statistic.
##  * p.value: the p-value of the test without bootstraping.
##  * p.value.boot: the p-value of the test with bootstraping.
##  * alternative: a string describing the type of alternative hypothesis.
##  * method: a string indicating what type of test was performed.
##  * data.name: a string giving the name(s) of the data.
##
mg_gof_ks.test2.boot <- function( x, y, n.boot = 1000, alternative = c("two.sided", "less", "greater"), exact = NULL, replace = FALSE )
{
	tol <- .Machine$double.eps*100;

	ks  <- ks.test( x, y, alt = alternative, ex = exact );

	n.x <- length(x);
	n.y <- length(y);
	cutp <- n.x;
	w <- c( x, y );
	n <- length(w);

	ks.pval <- 0; # estimated p-value
	#exact.boot <- ifelse( replace, FALSE, TRUE ); # with replacement we can get ties
	# Generates a new bunch of KS statistics under the null hypothesis
	# (i.e. x and y really come from the same distribution):
	# After the loop the p-value is estimated by:
	# <num of times KS-boot-stat >= KS-stat> / <num of boot>
	for ( k in 1:n.boot )
	{
		ix  <- sample( 1:n, n, rep = replace );

		xx <- w[ ix[1:cutp] ];
		yy <- w[ ix[(cutp+1):n] ];

		#ks.boot <- ks.test( xx, yy, alt = alternative, ex = exact.boot );
		ks.boot <- ks.test( xx, yy, alt = alternative, ex = exact );

		#if ( ks.boot$statistic >= (ks$statistic - tol) )
		if ( (ks$statistic - ks.boot$statistic) <= tol )
		{
			ks.pval  <- ks.pval + 1;
		}
	}
	ks.pval <- ks.pval/n.boot;

	return(
		list(
			statistic = ks$statistic,
			p.value = ks$p.value,
			p.value.boot = ks.pval,
			method = ks$method,
			alternative = ks$alternative,
			data.name = ks$data.name
		)
	);
}

#FIXME: incomplete
# Returns Pr( D_n <= d ), where D_n is the KS distribution.
mg_gof_ks.___cdf <- function( n, d )
{
	# taken from R source code "ks.test"
	pkstwo <- function(x, tol = 1e-06)
	{
		if (is.numeric(x)) 
		{
			x <- as.vector(x)
		} else {
			stop("argument 'x' must be numeric")
		}
		p <- rep(0, length(x))
		p[is.na(x)] <- NA
		IND <- which(!is.na(x) & (x > 0))
		if (length(IND) > 0) {
			p[IND] <- .C("pkstwo", as.integer(length(x[IND])), 
			p = as.double(x[IND]), as.double(tol), PACKAGE = "stats")$p
		}
		return(p)
	}
	return( pkstwo(sqrt(n) * d) );
}
