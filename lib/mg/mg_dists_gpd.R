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

## MG_DISTS_GPD
##
## SUMMARY
##  Generalized Pareto distribution
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

## MG_DISTS_DGPD
##
mg_dists_dgpd <- function( x, loc = 0, scale = 1, shape = 0, log = FALSE )
{
	if ( scale <= 0 )
	{
		return( rep( NaN, length(x) ) );
	}

	d <- (x - loc) / scale; # standardize
	good <- ( d >= 0 & ( (shape >= 0) | (shape*d > -1) ) ) | is.na(d); # x >= loc && (shape >=0 || loc <= x <= -(scale/shape) + loc)
	#if ( abs(shape) < .Machine$double.eps )
	if ( shape == 0 )
	{
		# Shape == 0

		d[good] <- - log( scale ) - d[good]; # log( (1/scale) * exp( - (x-loc)/shape ) )
	} else {
		# Shape != 0

		d[good] <- - log( scale ) - (1 + 1/shape) * log1p( shape * d[good] ); # log( (1/scale) * (1 + shape*((x-loc)/scale)^(-1-1/shape) )
	}
	d[!good] <- 0;

	if ( !log )
	{
		d[good] <- exp( d[good] );
	}

	return( d );
}

## MG_DISTS_PGPD
##
mg_dists_pgpd <- function( q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE )
{
	if ( scale <= 0 )
	{
		return( rep( NaN, length(q) ) );
	}

	p <- pmax( q - loc, 0 )/ scale; # standardize (for q < mu, force a zero value)

	if ( shape == 0 )
	{
		# Shape == 0

		p <- 1 - exp(-p); # 1 - exp( - (x-loc)/scale )
	} else {
		# Shape != 0

		p <- shape * p;
		good <- (shape >= 0) | (p > -1) | is.na(p); # x >= loc && (shape >=0 || loc <= q <= -(scale/shape) + loc)
		p[good] <- 1 - (1 + p[good])^(-1/shape); # 1 - ( 1 + shape*(q-loc)/scale )^(-1/shape)
		p[!good] <- 1; # (for (shape < 0 AND shape*(q-loc/scale) <= -1), force a one value)
	}

	if ( !lower.tail )
	{
		p <- 1 - p;
	}
	if ( log.p )
	{
		p <- log(p);
	}

	return( p );
}

mg_dists_qgpd <- function( p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE )
{
	if ( scale <= 0 )
	{
		return( rep( NaN, length(p) ) );
	}

	if ( log.p )
	{
		p <- log(p);
	}

	q <- rep( NaN, length(p) );

	good <- ( p >= 0 & p <= 1 ) | is.na(p);

	if ( shape == 0 )
	{
		# Shape == 0

		if ( lower.tail )
		{
			q[good] <- loc - scale * log1p( -p[good] );
		} else {
			q[good] <- loc - scale * log( p[good] );
		}
	} else {
		# Shape != 0

		if ( lower.tail )
		{
			q[good] <- loc + (scale / shape) * (1 - p[good])^shape - 1;
		} else {
			q[good] <- loc + (scale / shape) * ( p[good] )^shape - 1;
		}
	}

	return( q );
}

mg_dists_rgpd <- function( n, loc = 0, scale = 1, shape = 0 )
{
	if ( scale <= 0 )
	{
		return( rep( NaN, n ) );
	}

	return( mg_dists_qgpd( runif(n), loc = loc, scale = scale, shape = shape ) );
}

mg_dists_gpd.fit.mle <- function( x, loc = 0 )
{
	n.x <- length(x);
	x <- sort(x);
	x.range <- range(x);

	if ( n.x == 0 || !is.finite(x.range) )
	{
		return (
			list(
				loc = NaN,
				scale = NaN,
				shape = NaN
			)
		);
	} else if ( diff(x.range) < .Machine$double.xmin ) {
		# All observations are equal: try to return something reasonable
		if ( x.range[2] <= sqrt( .Machine$double.xmax ) )
		{
			shape <- NaN;
			scale <- 0;
		} else {
			shape <- -Inf;
			scale <- Inf;
		}
		return (
			list(
				loc = NaN,
				scale = scale,
				shape = shape
			)
		);
	}

	# Initial guess with MoM
	x.mean <- mean(x);
	x.var <- var(x);
	shape0 <- -0.5 * (x.mean^2 / x.var - 1);
	scale0 <- 0.5 * x.mean * (x.mean^2 / x.var + 1);
	if ( shape0 < 0 && (x.range[2] >= scale0/shape0 ) )
	{
		# MoM failed, start with an exponential fit
		shape0 <- 0;
		scale0 = x.mean;
	}
	params0 <- c( loc = loc, scale = scale0, shape = shape0 );

	negllik <- function( p, data ) {
		return( mg_dists_gpd.negllik( loc = p[1], scale = p[2], shape = p[3], data = data ) );
	};

	# Minimize the negative log-likelihood with respect to shape and ln(scale)
	res <- try( optim( params0, negllik, data = x ) );
	if ( is( res, "try-error" ) )
	{
		stop( "[MG::DISTS::GPD::FIT::MLE] Error while minimizing negative log-likelihood: ", geterrmessage() );
	}
	if ( res$convergence != 0 )
	{
		warning( "[MG::DISTS::GPD::FIT::MLE] Minimum negative log-likelihood estimation did not converge (exit code: ", res$convergence, ")" );
	}

	params.hat <- res$par;
	params.hat[2] <- params.hat[2];

	return(
		list(
			loc = params.hat[1],
			scale = params.hat[2],
			shape = params.hat[3]
		)
	);
}

## MG_DISTS_GPD.NEGLLIK( loc, scale, shape, x )
##
## SEE ALSO
##  * Computing Maximum Likelihood Estimates for the Generalized Pareto Distribution
##    Scott D. Grimshaw
##    Technometrics, Vol. 35, No. 2 (May, 1993), pp. 185-191
##    doi:10.2307/1269663
mg_dists_gpd.negllik <- function( loc, scale, shape, data )
{
	ln.scale <- log( scale ); # log-scale parameter
	n <- length(data);

	data[ data < loc | (shape == 0 & data > (loc - scale/shape)) ] <- 0;

	res <- Inf;
	if ( min( data ) >= loc )
	{
		# The support is x >= loc

		if ( shape == 0 )
		{
			# Shape == 0

			res <- n*ln.scale + sum( data - loc )/scale;
		} else {
			# Shape != 0

			if ( shape > 0 || max(data) <= (loc - scale/shape) )
			{
				# When shape < 0 the support is loc <= x <= (loc - scale/shape)
				# When shape >= the support is x >= loc

				res <- n*ln.scale + (1 + 1/shape) * sum( log1p( shape * (data - loc) / scale ) );
			}
		}
	}

	return( res );
}
