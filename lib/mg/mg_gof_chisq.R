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

## MG_GOF_CHISQ
##
## SUMMARY
##  Chi-Square Goodness-of-Fit tests.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

# TODO:
# * Chi-Square test for 2-sample (probabily we can use chisq.test)

source( "lib/mg/mg_dists.R" );

### MG_GOF_CHISQ.TEST
###
### SUMMARY
###  Chi-Square Goodness-of-Fit non-parametric test.
###
#mg_gof_chisq.test <- function( x, dist, ... )
#{
#	return( mg_gof_chisq.test.nonpar( x, dists, ... ) );
#}

## MG_GOF_CHISQ.TEST.PAR
##
## SUMMARY
##  Chi-Square Goodness-of-Fit parametric test.
##
## PARAMS
##  * x: data sample.
##  * dist: probability distribution name (@see MG_CONSTS).
##  * ...: parameters for @c dist.
##  * n.params: number of @c dist parameters estimated from sample @c x.
##  * e.freqs.min: minimun number of counts for each bin; if a bin count is less
##    than this number, bin pooling (merging) will be done.
##    Default value is 5. Use 0 to prevent pooling.
##  * bin.method: method for binning the distribution ("moore" or "sturges").
##
## RETURN
##  A list with the following field:
##  * statistic: the Chi-Square test statistic.
##  * parameter: the number of degrees of freedom for Chi-Square distribution.
##  * p.value: the p-value for the test.
##  * observed: the observed frequencies.
##  * expected: the expected frequencies.
##
## DESCRIPTION
##  This is a parametric version of Chi-Square test. Use this when the
##  theoretical distribution @c dist has some (eventually all) parameters
##  estimated from the sample @c x against which you want to perform the test.
##
mg_gof_chisq.test.par <- function( x, dist, ..., n.params, e.freqs.min = 5, bin.method = c("moore", "sturges") )
{
	return( mg_gof_chisq._test( x, dist, ..., n.params, e.freqs.min = e.freqs.min, bin.method = bin.method ) );
}

## MG_GOF_CHISQ.TEST.NONPAR
##
## SUMMARY
##  Chi-Square Goodness-of-Fit non-parametric test.
##
## PARAMS
##  * x: data sample.
##  * dist: probability distribution name (@see MG_CONSTS).
##  * ...: parameters for @c dist.
##  * e.freqs.min: minimun number of counts for each bin; if a bin count is less
##    than this number, bin pooling (merging) will be done.
##    Default value is 5. Use 0 to prevent pooling.
##  * bin.method: method for binning the distribution ("moore" or "sturges").
##
## RETURN
##  A list with the following field:
##  * statistic: the Chi-Square test statistic.
##  * parameter: the number of degrees of freedom for Chi-Square distribution.
##  * p.value: the p-value for the test.
##  * observed: the observed frequencies.
##  * expected: the expected frequencies.
##
## DESCRIPTION
##  This is a non-parametric version of Chi-Square test. Use this when no one
##  parameter of the theoretical distribution @c dist has been estimated from
##  the sample @c x.
##
mg_gof_chisq.test.nonpar <- function( x, dist, ..., e.freqs.min = 5, bin.method = c("moore", "sturges") )
{
	return( mg_gof_chisq._test( x, dist, ..., e.freqs.min = e.freqs.min, bin.method = bin.method ) );
}

## MG_GOF_CHISQ.TEST
##
## SUMMARY
##  Chi-Square Goodness-of-Fit test.
##
## PARAMS
##  * x: data sample.
##  * dist: probability distribution name (@see MG_CONSTS).
##  * ...: parameters for @c dist.
##  * n.params: number of @c dist parameters estimated from sample @c x.
##  * e.freqs.min: minimun number of counts for each bin; if a bin count is less
##    than this number, bin pooling (merging) will be done.
##    Default value is 5. Use 0 to prevent pooling.
##  * bin.method: method for binning the distribution ("moore" or "sturges").
##
## RETURN
##  A list with the following field:
##  * statistic: the Chi-Square test statistic.
##  * parameter: the number of degrees of freedom for Chi-Square distribution.
##  * p.value: the p-value for the test.
##  * observed: the observed frequencies.
##  * expected: the expected frequencies.
##
## DESCRIPTION
##  This is an internal function used to perform either parametric and
##  non-parametric Chi-Square test.
##
mg_gof_chisq._test <- function( x, dist, ..., n.params = 0, e.freqs.min = 5, bin.method = c("moore", "sturges") )
{
	bin.method <- match.arg( bin.method );
	x.len <- length(x);
	n.bins <- 0;
	bin.breaks <- c();

	if ( mg_dists_isDiscrete( dist ) )
	{
		## >>> TODO: this part need to tested!!! <<<

		# Creates bins
		xx <- sort(unique(x));
 		#x.breaks <- c(floor(xx[1])-1, xx, ceiling(xx[length(xx)])+1);
 		bin.breaks <- c( xx, ceiling(xx[length(xx)])+1);
		o.hist <- hist( x, breaks = bin.breaks, right = FALSE, plot = FALSE );

		# Computes observed frequencies
		o.freqs <- o.hist$counts;

		# Computes expected frequencies
		#e.freqs <- diff( c( 0, mg_dists_cdf( dist, bin.breaks[1:length(o.freqs)], ... ) ) ) * x.len;
		e.freqs <- try( mg_dists_cdf( dist, bin.breaks[2:(length(o.freqs)-1)], ... ) );
		if ( is( e.freqs, "try-error" ) )
		{
			stop( "[MG::GOF::CHISQ::_TEST] Error while evaluating CDF for distribution '", dist, "': ", geterrmessage() );
		}
		e.freqs <- diff( c( 0, e.freqs, 1 ) ) * sum(o.freqs);

		n.bins <- length( o.freqs );
	} else {
##@{ Deprecated [2007-05-25]
#		# Divide the distribution support (i.e. range) in bins.
#
#		x.len <- length(x);
#		#rvg <- get( mg_dists_rvgName( dist ), mode = "function" );
#		#y <- rvg( x.len, ... );
#		y <- mg_dists_rvg( dist, x.len, ... );
#		hy <- hist( y, plot = FALSE );
#		x.cut <- cut( x, breaks = length(hy$density) )
#		#hx <- hist( x, breaks = length(hy$breaks), plot = FALSE );
#		#res$chisq <- chisq.test( h$mid, h$density );
#		return( chisq.test( (table(x.cut)/x.len)[], hy$density ) );
##@} Deprecated [2007-05-25]

		#x.len <- length(x);

		# Creates bins
		n.bins <- 0;
		if ( bin.method == "sturges" )
		{
			n.bins <- ceiling( 1 + log2(x.len) ); #see http://www.mathwave.com/articles/goodness_of_fit.html
		} else if ( bin.method == "moore" )
		{
			n.bins <- ceiling( 2*x.len^(2/5) ); # See Mann-Wald (1942), Moore (1986) and Del Barrio (2000)
		} else {
			stop( "[MG::GOF::CHISQ::_TEST] Unknown binning method '", bin.method, "'" );
		}
		bin.area <- 1/n.bins;
		#p <- bin.area * c(0:n.bins);
		p <- bin.area * (1:n.bins - 0.5);
		bin.breaks <- try( mg_dists_invcdf( dist, p, ... ) );
		if ( is( bin.breaks, "try-error" ) )
		{
			stop( "[MG::GOF::CHISQ::_TEST] Error while evaluating quantile function for distribution '", dist, "': ", geterrmessage() );
		}
		#bin.breaks[is.infinite(bin.breaks)] <- max(x); # handle +
		x.range <- range(x);
		bin.breaks[which(bin.breaks == -Inf)] <- x.range[1]; # handle -Inf
		bin.breaks[which(bin.breaks == Inf)] <- x.range[2]; # handle +Inf
		#bin.breaks <- mg_dists_chisq._createbins( x, n.bins, byprob = TRUE );

		# Computes observed frequencies for each bin
		##intervals <- findIntervals( x, bin.breaks )
		##o.freqs <- hist( intervals, breaks = c( (intervals[1] - 0.5):(intervals[length(intervals)] + 0.5) ) );
		#o.freqs <- hist( x, breaks = c( -Inf, bin.breaks[2:(length(bin.breaks)-1)], +Inf ), plot = FALSE )$counts; # No! hist don't handle Inf values
		#o.freqs <- hist( x, breaks = c( -.Machine$double.xmax, bin.breaks[2:(length(bin.breaks)-1)], +.Machine$double.xmax ), plot = FALSE )$counts;
		breaks.range <- range(bin.breaks);
		breaks.len <- length(bin.breaks);
		o.freqs <- hist(
			x,
			breaks = c( min( x.range[1], breaks.range[1] ), bin.breaks[2:(breaks.len-1)], max( x.range[2], breaks.range[2] ) ),
			plot = FALSE,
			right = FALSE
		)$counts;

		# Computes expected frequencies: e.freqs[i] <- n * (F(breaks[2]) - F(breaks[1]))
		#e.freqs.old <- rep( mg_dists_cdf( dist, bin.breaks[2], ... ) - mg_dists_cdf( dist, bin.breaks[1], ... ), n.bins-1 ) * x.len;
		e.freqs <- try( mg_dists_cdf( dist, bin.breaks[2:(breaks.len-1)], ... ) );
		if ( is( e.freqs, "try-error" ) )
		{
			stop( "[MG::GOF::CHISQ::_TEST] Error while evaluating CDF for distribution '", dist, "': ", geterrmessage() );
		}
		e.freqs <- sum( o.freqs ) * diff( c( 0, e.freqs, 1 ) );
	}

	o.freqs[ is.na(o.freqs) ] <- 0;
	e.freqs[ is.na(e.freqs) ] <- 0;

	# Checks for pooling
	if ( any(e.freqs) < e.freqs.min )
	{
		# Performs pooling

		pooled <- mg_gof_chisq._poolbins( o.freqs, e.freqs, bin.breaks, e.freqs.min )
		o.freqs <- pooled$o.freqs;
		e.freqs <- pooled$e.freqs;
		bin.breaks <- pooled$breaks;
		n.bins <- length( o.freqs );
	}

	# Computes the test statistics
	chisq.stat <- sum( (o.freqs - e.freqs)^2 / e.freqs );
	#chisq.stat <- sum( o.freqs^2/e.freqs ) - x.len; # a more efficient formula
	# Sets the degrees of freedom
	chisq.df <- n.bins - 1 - n.params;

	# Computes the p-value
	if ( chisq.df > 0 )
	{
		chisq.pval <- 1 - pchisq( chisq.stat, chisq.df );
	} else {
		chisq.df <- 0;
		chisq.pval <- NaN;
	}

	return( 
		list(
			statistic = chisq.stat,
			parameter = chisq.df,
			p.value = chisq.pval,
			observed = o.freqs,
			expected = e.freqs,
			residuals = (o.freqs - e.freqs) / sqrt(e.freqs)
			#data.namme = ...
			#method = ...
		)
	);
}

#mg_gof_chisq._createbins <- function( x, nbins, byprob = TRUE )
#{
#	bin.breaks <- c();
#
#	if ( byprob )
#	{
#		x.range <- range(x);
#		if ( is.null(x) || length(x) == 0 )
#		{
#			x.range[1] <- 0;
#			x.range[2] <- 1;
#		}
#		if ( x.range[1] == x.range[2] )
#		{
#			x.range[1] <- x.range[1] - floor( nbins/2 ) - 0.5;
#			x.range[2] <- x.range[2] - ceiling( nbins/2 ) - 0.5;
#		}
#
#		# Creates edges and centers
#		bin.width <- diff(x.range) / nbins;
#		edges <- range[1] + bin.width * (0:nbins);
#		edges.len <- length(edges);
#		edges[edges.len] <- x.range[2];
#		mids <- edges[ 1:(length(edges)-1) ] + bin.width/2;
#
#		# Update bin widths for internal bins
#		nbins <- edges.len - 1;
#
#		# Shift bins so the internal is ( ] instead of [ ).
#		edges <- edges + .Machine$double.eps;
#
#		# Map each jump location to a bin number. -Inf accounts for the above
#		# shift, +Inf keeps things out of histc's degenerate rightmost bin.
#		bin.breaks <- c( -Inf, edges[2:(edges.len-1)] Inf );
#	}
#
#	return( bin.breaks );
#}

## MG_GOF_CHISQ._POOLBINS
## 
## SUMMARY
##  Check that expected bin frequencies are not too small.
##
## DESCRIPTION
##  Pool the smallest bin each time, working from the end, but
##  avoid pooling everything into one bin.  We will never pool bins
##  except at either edge (no two internal bins will get pooled together).
##
mg_gof_chisq._poolbins <- function( o.freqs, e.freqs, breaks, e.freqs.min )
{
	# taken from MATLAB R2006b

	i <- 1;
	j <- length( e.freqs );
	# Enters the loop only if two adjacent bins (one of which has freqs < freqs.min) are found
	while ( i < (j-1) && ( e.freqs[i] < e.freqs.min || e.freqs[i+1] < e.freqs.min || e.freqs[j] < e.freqs.min || e.freqs[j-1] < e.freqs.min ) )
	{
		# Found two adjacent bins
		if ( e.freqs[i] < e.freqs[j] )
		{
			e.freqs[i+1] <- e.freqs[i+1] + e.freqs[i];
			o.freqs[i+1] <- o.freqs[i+1] + o.freqs[i];
			i <- i+1;
		} else {
			e.freqs[j-1] <- e.freqs[j-1] + e.freqs[j];
			o.freqs[j-1] <- o.freqs[j-1] + o.freqs[j];
			j <- j-1;
		}
	}

	# Retain only the pooled bins
	e.freqs <- e.freqs[i:j];
	o.freqs <- o.freqs[i:j];
	if ( j < (length(breaks) - 1) )
	{
		# removes (j+1):(length(breaks)+1) breaks
		#breaks <- c( breaks[ 1:j ], breaks[length(breaks)] );
		breaks <- breaks[ -( (j+1):(length(breaks)-1) ) ];
	}
	if ( i > 1 )
	{
		# removes 2:i breaks
		#breaks <- c( breaks[1], breaks[(i+1):length(breaks)] );
		breaks <- breaks[ -( 2:i ) ];
	}

	# Warn if some remaining bins have expected counts too low
	if ( any( e.freqs < e.freqs.min ) )
	{
		warning( "[MG::GOF::CHISQ] After pooling, some bins still have low expected counts.\n", "The chi-square approximation may be inaccurate." );
	}

	return(
		list(
			o.freqs = o.freqs,
			e.freqs = e.freqs,
			breaks = breaks
		)
	);
}
