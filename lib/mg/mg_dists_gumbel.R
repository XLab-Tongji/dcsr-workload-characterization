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

## Gumbel (Type 1 Extreme Value) distribution
##
## SUMMARY
##  Gumbel or Type I Extreme Value distribution.
##
## DESCRIPTION
##  Gumbel or Type I Extreme Value distribution:
##   f(x)=\frac{\exp{-\frac{x-\mu}{\sigma}}\exp(-\exp{-\frac{x-\mu}{\sigma}})}{\sigma}
##   F(x;\mu,\sigma)=\exp{-\exp{-(x-\mu)/\sigma}} x\in\mathbb{R}
##
##  Distribution Parameters:
##  * mu: location parameter (range (-Inf,+Inf)
##  * sigma: scale parameter (range (0,+Inf)
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

mg_dists_dgumbel <- function( x, location = 0, scale = 1, log = FALSE )
{
#	mg_dists_gumbel_haveevd <- require( "evd" );
#
#	if ( mg_dists_gumbel_haveevd )
#	{
#		return( evd::dgumbel( x, loc = location, scale = scale ) );
#	}

	if ( any( scale <= 0 ) )
	{
		stop( "[MG::Dists::Gumbel::DGumbel] Scale parameter out-of-range (<=0)." );
		#return( seq( NaN, length(x) ) );
	}

	z <- ( x - location ) / scale;
	d <- -log(scale) - z - exp( - z );

	if ( !log )
	{
		d <- exp( d );
	}

	return( d );
}

mg_dists_pgumbel <- function( q, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE )
{
#	mg_dists_gumbel_haveevd <- require( "evd" );
#
#	if ( mg_dists_gumbel_haveevd )
#	{
#		p <- evd::pgumbel( q, loc = location, scale = scale, lower.tail );
#	} else {
		if ( any( scale <= 0 ) )
		{
			stop( "[MG::Dists::Gumbel::PGumbel] Scale parameter out-of-range (<=0)." );
		}

		z <- ( q - location ) / scale;
		p <- exp( - exp( -z ) );

		if ( !lower.tail )
		{
			p <- 1 - p;
		}
#	}

	if ( log.p )
	{
		p <- log( p );
	}

	return( p );
}

mg_dists_qgumbel <- function( p, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE )
{
#	mg_dists_gumbel_haveevd <- require( "evd" );
#
#	if ( mg_dists_gumbel_haveevd )
#	{
#		if ( log.p )
#		{
#			p <- exp( p );
#		}
#		return( evd::qgumbel( p, loc = location, scale = scale, lower.tail ) );
#	}

	if ( any( scale <= 0 ) )
	{
		stop( "[MG::Dists::Gumbel::QGumbel] Scale parameter out-of-range (<=0)." );
	}

	if ( log.p )
	{
		p <- exp( p );
	}
	if ( !lower.tail )
	{
		p <- 1 - p;
	}

	return( location - scale * log( -log( p ) ) );
}

mg_dists_rgumbel <- function( n, location = 0, scale = 1 )
{
#	mg_dists_gumbel_haveevd <- require( "evd" );
#
#	if ( mg_dists_gumbel_haveevd )
#	{
#		return( evd::rgumbel( n, loc = location, scale = scale ) );
#	}

	if ( any( scale <= 0 ) )
	{
		stop( "[MG::Dists::Gumbel::RGumbel] Scale parameter out-of-range (<=0)." );
	}

	return( location - scale * log( - log( runif( n ) ) ) ); # Alternative #1
	#return( location - scale * log( rexp( n ) ) ); # Alternative #2
}

## Fits a Gumbel distribution with method of moments
mg_dists_gumbel.fit.mm <- function( x )
{
	sigmahat <- sqrt(6) * sd(x) / pi;
	muhat <- mean(x) - 0.57722 * sigma0;

	return( list( location = muhat, scale = sigmahat ) );
}

## Fits a Gumbel distribution with Maximum Likelihood method
mg_dists_gumbel.fit.mle <- function( x )
{
	n <- length(x);

	# Use method-of-moment for estimating initial values
	sigma0 <- sqrt(6) * sd(x) / pi;
	mu0 <- mean(x) - 0.57722 * sigma0;

	parms0 <- c( sigma0, mu0 );

	# The negative log-likelihood function
	negloglik <- function( parmshat, x )
	{
		# parmshat == c( sigma, mu )

		if ( parmshat[1] <= 0 )
		{
			res <- 1e+06;
		} else {
			y <- (x - parmshat[2]) / parmshat[1];

			term1 <- length(x) * log( parmshat[1] ); 
			#term2 <- sum( x ) - length(x)*parmshat[2]/parmshat[1];
			term2 <- sum( y );
			term3 <- sum( exp( - y ) );
			res <- term1 + term2 + term3;
		}
		return( res );
	};

	# Solve log-Likelihood(parmshat) = 0
	fit <- optim( parms0, negloglik, x = x );
	if ( fit$convergence )
	{
		warning(" [MG::Dists::Gumbel::Fit::MLE] Optimization may not have succeeded." );
	}
	return( list( location = fit$par[2], scale = fit$par[1] ) );
}

# This is the implementation used in MATLAB (see evfit).
# The R implementation is consistent with MATLAB one (i.e. same
# result); however both R and MATLAB implementation seems not to give
# the correct result (e.g. compares the results with ones obtained with
# the above mg_dists_gumbel.fit.mle function).
mg_dists_gumbel.fit.mle.DontWork <- function( x )
{
	n <- length( x );
	rangex <- diff( range( x ) );
	maxx <- max( x );
	xtol <- .Machine$double.eps^0.25;

	if ( rangex < .Machine$double.xmin )
	{
		# The likelihood surface for constant data has its maximum at the
		# boundary sigma=0, Return something reasonable anyway.
		return( list( location = x[1], scale = 0 ) );
	}

	lkeqn <- function(sigma, x, xbarWgtUnc)
	{
		# Likelihood equation for the extreme value scale parameter.
		# Assumes that sigma is strictly positive and finite, and x has
		# been transformed to be non-positive and have a reasonable
		# range.
		#
		# as sigma->0, v->(sigma+xbarWgtUnc-0)->xbarWgtUnc < 0
		# as sigma->Inf, v ->(sigma+xbarWgtUnc-xbar)->Inf > 0
		#
		w <- exp(x / sigma);
		v <- sigma + xbarWgtUnc - sum(x * w) / sum(w);

		return( v );
	};

	# Shift x to max(x) == 0, min(x) = -1 to make likelihood eqn more stable.
	x0 <- (x - maxx) / rangex;

	# First, get a rough estimate for the scale parameter sigma as a starting value.
	# (use MM)
	sigmahat <- (sqrt(6) * sd(x0)) / pi;
	#wgtmeanUnc <- sum(x0) / n;
	wgtmeanUnc <- mean(x0);

	# Bracket the root of the scale parameter likelihood eqn ...
	if ( lkeqn(sigmahat, x0, wgtmeanUnc) > 0 )
	{
		upper <- sigmahat;
		lower <- .5 * upper;
		while ( lkeqn(lower, x0, wgtmeanUnc) > 0 )
		{
			upper <- lower;
			lower <- .5 * upper;
			if ( lower < .Machine$double.xmin ) #  underflow, no positive root
			{
				stop("[MG::FIT::MLEGUMBEL] No maximum likelihood solution found.");
			}
		}
	} else {
		lower <- sigmahat;
		upper <- 2 * lower;
		while ( lkeqn(upper, x0, wgtmeanUnc) < 0 )
		{
			lower <- upper;
			upper <- 2 * lower;
			if ( upper > .Machine$double.xmax ) # overflow, no finite root
			{
				stop("[MG::FIT::MLEGUMBEL] No maximum likelihood solution found.");
			}
		}
	}
	bnds <- c(lower, upper);

	# ... then find the root of the likelihood eqn.  That is the MLE for sigma,
	# and the MLE for mu has an explicit sol'n. 
	r <- uniroot( lkeqn, bnds, tol = xtol, x = x0, xbarWgtUnc = wgtmeanUnc );
	sigmahat <- r$root;
	if ( abs(r$f.root) > xtol )
	{
		warning("[MG::FIT::MLEGUMBEL] The likelihood equation may be ill-conditioned.");
	}
	#muhat <- sigmahat * log( sum( exp(x0/sigmahat) ) / n );
	muhat <- sigmahat * log( mean( exp(x0/sigmahat) ) );

	# Those were parameter estimates for the shifted, scaled data, now
	# transform the parameters back to the original location and scale.
	return( list( location = (rangex*muhat)+maxx, scale = rangex*sigmahat ) );
}
