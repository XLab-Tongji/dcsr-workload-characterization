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

## MG_MATH_FUNCS
##
## SUMMARY
##  A collection of special mathematics functions.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

##
## MG_MATH_LOG2
##
## SUMMARY
##  Base 2 logarithm and dissect floating-point numbers into exponent and
##  mantissa.
##
## DESCRIPTION
##  Returns arrays @c l, @c f and @c e.
##
##  Return value @c l is an array of real values representing the base 2 log of
##  @c x.
##  Return value @c f is an array of real values, usually in the range
##  @c 0.5 <= abs(f) < 1. For real @c x, @c f satisfies the equation:
##  @c x = f.*2.^e.
##  Return value @c e is an array of integers that, for real @c x, satisfy the
##  equation: @c x = f.*2.^e.
##
## REMARKS
##  This function corresponds to the ANSI C function frexp() and the IEEE
##  floating-point standard function logb(). Any zeros in X produce F = 0 and
##  E = 0.
##
mg_math_log2 <- function( x )
{
	ll <- log2(x); # The base-2 log values

	x <- Re(x); # only deals with the real parts ...

	# Note: 0 entries are replaced by 1, by multiplication with the sign.

	ff <- abs(x) + (x == 0);
	ee <- floor( log(ff)/log(2) + 1 ) * (x != 0);
	ff <- sign(x) * ff / ( 2^ee );

	return( list( l = ll, f = ff, e = ee ) );
}
