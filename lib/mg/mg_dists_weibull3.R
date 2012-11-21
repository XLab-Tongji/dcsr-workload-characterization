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

mg_dists_dweibull3 <- function( x, shape, scale = 1, loc = 0, log = FALSE )
{
	# Looks for out-of-range values
	good <- x >= shape;

	# Out-of-range values take a zero
	d <- rep.int( 0, length(x) );

	# Use a two-parameter Weibull with a change of variable Y=(X-loc)
	d[good] <- try(
		dweibull(
			x[good] - loc,
			shape = shape,
			scale = scale,
			log = log
		)
	);
	if ( any( is( d[good], "try-error" ) ) )
	{
		stop( "[MG::DISTS::DWEIBULL3] Error while calculating densities: ", geterrmessage(), "." );
	} 

	return( d );
}

mg_dists_pweibull3 <- function( q, shape, scale = 1, loc = 0, lower.tail = TRUE, log.p = FALSE )
{
	# Looks for out-of-range values
	good <- q >= shape;

	# Out-of-range values take a zero
	p <- rep.int( 0, length(q) );

	# Use a two-parameter Weibull with a change of variable Y=(X-loc)
	p[good] <- try(
		pweibull(
			q[good] - loc,
			shape = shape,
			scale = scale,
			lower.tail = lower.tail,
			log.p = log.p
		)
	);
	if ( any( is( p[good], "try-error" ) ) )
	{
		stop( "[MG::DISTS::PWEIBULL3] Error while calculating probabilities: ", geterrmessage(), "." );
	} 

	return( p );
}

mg_dists_qweibull3 <- function( p, shape, scale = 1, loc = 0, lower.tail = TRUE, log.p = FALSE )
{
	q <- rep.int( NaN, length(q) );

	# Looks for out-of-range values
	good <- p >= 0 & p <= 1;

	# Use a two-parameter Weibull with a change of variable Y=(X-loc)
	q[good] <- try(
		qweibull(
			p[good],
			shape = shape,
			scale = scale,
			lower.tail = lower.tail,
			log.p = log.p
		)
		+ loc
	);
	if ( any( is( q[good], "try-error" ) ) )
	{
		stop( "[MG::DISTS::QWEIBULL3] Error while calculating quantiles: ", geterrmessage(), "." );
	} 

	# Out-of-range values take a infinity
	q[p < 0] <- -Inf;
	q[p > 1] <- Inf;

	return( q );
}

mg_dists_rweibull3 <- function( n, shape, scale = 1, loc = 0 )
{
	# Use a two-parameter Weibull
	x <- try(
		rweibull(
			n,
			shape = shape,
			scale = scale
		)
		+ loc
	);
	if ( any( is( x, "try-error" ) ) )
	{
		stop( "[MG::DISTS::RWEIBULL3] Error while generating random numbers: ", geterrmessage(), "." );
	} 

	return( x );
}

mg_dists_weibull3.fit.mle <- function( x )
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
                        scale = as.numeric( x$estimate["scale"] ),
			loc = min(x) # FIXME!!!!!!!!
                )
        );
}

mg_dists_weibull.fit.mle.test <- function( x )
{
        n <- length(x);

	m <- mean( log(x) );
	var <- var( log(x) );
	loc0 <- min(x);
	shape0 <- 1.2 / sqrt(v);
	scale0 <- exp( m + 0.572/shape0 );

	param.hat.ev <- evfit( log(x) );
	param.hat <- c( exp( param.hat.ev[1] ), 1/param.hat.ev[2] );
}
