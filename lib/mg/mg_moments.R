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

## MG_MOMENTS
##
## SUMMARY
##  A collection of moments utilities.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

## Returns the k-th sample moments; if "central" is TRUE returns the k-th
## central sample moment.
mg_moments_moment <- function( x, k, central = FALSE )
{
	if ( central )
	{
		return( mean( ( x - mean( x ) )^k ) );
	}

	return( mean( x^k ) );
}

## Returns the skewness (3rd standardized - or normalized - moment)
##   gamma_1 = mu_3 / mu_2^(3/2)
## where mu_3 is the 3rd central moment about the mean and mu_2=sigma is the
## standard deviation (2nd central moment about the mean).
mg_moments_skewness <- function( x )
{
	#return( mg_moments_moment( x, 3, TRUE ) / ( sqrt( mg_moments_moment( x, 2, TRUE ) )^3 ) );
	return( mg_moments_moment( x, 3, TRUE ) / ( sd(x)^3 ) );
}

## Returns the Quartile coefficient of skewness.
## Positive skewness indicates that the median is closer to the
## lower quartile than the upper quartile, negative that the median is closer to
## the upper quartile than the lower quartile.
## It's value is in [0,1].
mg_moments_skewness.quartile <- function( x )
{
	qq <- quantile( x, probs = c(0.25, 0.5, 0.75) );

	#qcsk <- (x25+x75-2*mean(x))/(x75-x25)
	return( (qq[1] - 2*qq[2] + qq[3]) / (qq[3] - qq[1]) );
}

## Returns the kurtosis coefficient.
## This is equivalent to the fourth central moment divided by the variance squared.
## Kurtosis > 3 indicates a distribution that is more peaked and has heavier tails
## than a normal distribution with the same variance.  Kurtosis < 3 indicates a
## distribution that is flatter or has heavier flanks than the normal.
mg_moments_kurtosis <- function( x )
{
	return( mg_moments_moment( x, 4, TRUE ) / var(x)^2 );
}

## Returns the kurtosis excess.
## Positive kurtosis excess indicates a distribution that is more peaked and has
## heavier tails than a normal distribution with the same variance.  Negative
## kurtosis excess indicates a distribution that is flatter or has heavier
## flanks than the normal.
mg_moments_kurtosis.excess <- function( x )
{
	return( mg_moments_kurtosis( x ) - 3 );
}
