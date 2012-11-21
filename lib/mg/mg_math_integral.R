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

## MG_MATH_INTEGRAL
##
## Function for solving numerical integrals.
## 
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

## Calculates the R default integration method
mg_math_integral.quad <- function( f, lower, upper, subdivisions = 100, rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol, stop.on.error = TRUE, keep.xy = FALSE, aux = NULL, ... )
{
	return( integrate( f, lower, upper, subdivisions, rel.tol, abs.tol, stop.on.error, keep.xy, aux, ... ) );
}

##
## Calculates the numerical integral with the simple mid-point (or rectangle)
## formula:
##  \int_l^u f(x)dx ~= (u-l) * f( (l+u)/2 )
## where:
## * "f" is the function to be integrated
##
mg_math_integral.midpoint <- function( f, lower, upper, ... )
{
	f <- match.fun( f );
	ff <- function(x) { return( f( x, ... ) ); };

	return( (upper - lower) * ff( (lower + upper)/2 ) );
}

##
## Calculates the cumulative numerical integral with the composite mid-point (or
## rectangle) formula:
##  \int_l^u f(x)dx ~= H * \sum_{k = 1}^M f( \bar{x}_k )
## where:
## * "f" is the function to be integrated
## * "H" is the sub-interval size ("H = (u-l)/M", where "l" is the lower limit
##   and "u" is the upper limit)
## * "M" is the number of subdivisions (i.e. sub-intervals).
##
mg_math_integral.midpointc <- function( f, lower, upper, subdivisions = 100, ... )
{
	f <- match.fun( f );
	ff <- function(x) { return( f( x, ... ) ); };

	limit <- as.integer( subdivisions );
	H <- (upper - lower) / limit;

	return(
		H * sum(
			ff(
				seq(
					lower + H/2,
					upper - H/2,
					length = limit
				)
			)
		)
	);
}

##
## Calculates the numerical integral with the simple trapezoidal formula:
##  \int_l^u f(x)dx ~= (u-l)/2 * \[f(l)+f(u)\]
## where:
## * "f" is the function to be integrated
##
mg_math_integral.trapz <- function( f, lower, upper, subdivisions = 100, ... )
{
	f <- match.fun( f );
	ff <- function(x) { return( f( x, ... ) ); };

	return( (upper - lower)/2 * (ff(lower) + ff(upper)) );
}

##
## Calculates the cumulative numerical integral with the composite trapezoidal
## formula:
##  \int_l^u f(x)dx ~= H/2 * \sum_{k = 1}^M f(x_{k-1})+f(x_k)
## where:
## * "f" is the function to be integrated
## * "H" is the sub-interval size ("H = (u-l)/M", where "l" is the lower limit
##   and "u" is the upper limit)
## * "M" is the number of subdivisions (i.e. sub-intervals).
##
mg_math_integral.trapzc <- function( f, lower, upper, subdivisions = 100, ... )
{
	f <- match.fun( f );
	ff <- function(x) { return( f( x, ... ) ); };

	limit <- as.integer( subdivisions );
	H <- (upper - lower) / limit;

	xx <- seq( lower, upper, length = limit + 1 );
#	yy <- ff(xx);
#	yy[2:limit] <- 2 * yy[2:limit]
#
#	return( 0.5 * H * sum(yy) );

	return( 0.5 * H * (ff(xx[1]) + ff(xx[limit+1])) + H * sum( ff(xx[2:limit]) ) );
}

##
## Calculates the numerical integral with the 2-nodes simple gaussian quadrature
## formula:
##  \int_l^u f(x)dx ~= (u-l)/2 * \[f(\gamma_0)+f(\gamma_1)\]
## where:
## * "f" is the function to be integrated
## * "\gamma_0 = x_{k-1} + (1-1/sqrt(3))*H/2" is a Gauss node.
## * "\gamma_1 = x_{k-1} + (1+1/sqrt(3))*H/2" is a Gauss node.
##
mg_math_integral.gauss2 <- function( f, lower, upper, subdivisions = 100, ... )
{
	f <- match.fun( f );
	ff <- function(x) { return( f( x, ... ) ); };

	h <- (upper - lower) / (2 * sqrt(3));
	gg0 <- (lower + upper)/2 - h; # gamma 0
	gg1 <- (lower + upper)/2 + h; # gamma 1

	return( (upper - lower)/2 * (ff(gg0) + ff(gg1) ) );
}

##
## Calculates the cumulative numerical integral with the 2-nodes composite
## gaussian quadrature formula:
##  \int_l^u f(x)dx ~= H/2 * \sum_{k = 1}^M f(\gamma_{k-1})+f(\gamma_k)
## where:
## * "f" is the function to be integrated
## * "\gamma_{k-1} = x_{k-1} + (1-1/sqrt(3))*H/2" is a Gauss node.
## * "\gamma_k = x_{k-1} + (1+1/sqrt(3))*H/2" is a Gauss node.
## * "H" is the sub-interval size ("H = (u-l)/M", where "l" is the lower limit
##   and "u" is the upper limit)
## * "M" is the number of subdivisions (i.e. sub-intervals).
##
mg_math_integral.gauss2c <- function( f, lower, upper, subdivisions = 100, ... )
{
	f <- match.fun( f );
	ff <- function(x) { return( f( x, ... ) ); };

	limit <- as.integer( subdivisions );
	H <- (upper - lower) / limit;

	xx <- seq( lower, upper, length = limit + 1 );
	yy <- xx[1:limit] + H/2 * (1 - 1/sqrt(3));
	yy <- c( yy, yy + H/sqrt(3) );

	return( 0.5 * H * sum( ff(yy) ) );
}

##
## Calculates the numerical integral with the simple Simpson formula:
##  \int_l^u f(x)dx ~= (u-l)/6 * \[f(l)+4*f((l+u)/2)+f(u)\]
## where:
## * "f" is the function to be integrated.
##
mg_math_integral.simpson <- function( f, lower, upper, ... )
{
	f <- match.fun( f );
	ff <- function(x) { return( f( x, ... ) ); };

	return( (upper-lower)/6 * (ff(lower) + 4*ff((upper+lower)/2) + ff(upper) ) );
}

##
## Calculates the cumulative numerical integral with the composite Simpson
## formula:
##  \int_l^u f(x)dx ~= H/6 * \sum_{k = 1}^M \[f(x_{k-1})+4*f(\bar{x}_k)+f(\x_k)\]
## where:
## * "f" is the function to be integrated
## * "H" is the sub-interval size ("H = (u-l)/M", where "l" is the lower limit
##   and "u" is the upper limit)
## * "M" is the number of subdivisions (i.e. sub-intervals).
##
mg_math_integral.simpsonc <- function( f, lower, upper, subdivisions = 100, ... )
{
	f <- match.fun( f );
	ff <- function(x) { return( f( x, ... ) ); };

	limit <- as.integer( subdivisions );
	H <- (upper - lower) / limit;

	xx <- seq( lower, upper, length = limit + 1 );
	yy <- ff( xx );
	yy[2:limit] <- 2 * yy[2:limit];
	yy <- H * sum( yy ) / 6;
	xx <- seq( lower + H/2, upper - H/2, length = limit )

	return( yy + 2*H*sum( ff(xx) )/3 );
}
