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

## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

## MG_PLOT_LRD.AGGVAR
##
mg_plot_lrd.aggvar <- function( av, n.levels = NULL, min.points = 3, cutoffs = 10^c(0.7, 2.5), main = NULL, xlab = "m", ylab = "variance", col = "royalblue", plot.refline = TRUE, plot.fitline = FALSE, plot.cutlines = TRUE )
{
	if (
		!is( av, "mg_lrd_aggvar" )
	)
	{
		av <- as.vector(av);
		av <- mg_lrd_aggvar(
			av,
			n.levels = ifelse( !is.null(n.levels), n.levels, log10(length(av)) ),
			min.points = min.points,
			cutoffs = cutoffs,
			ret.m = TRUE,
			ret.vars = TRUE,
			ret.log.m = TRUE,
			ret.log.vars = TRUE
		);
	}
	if (
		( is.null(av$m) && is.null(av$vars) )
		&& ( is.null(av$log.m) && is.null(av$log.vars) )
	)
	{
		stop( "[MG::PLOT::LRD::AGGVAR] Cannot plot Aggregated Variances without (x,y) values." );
	}

	if ( is.null( av$log.m ) )
	{
		av$log.m <- log10( av$m );
		av$log.vars <- log10( av$vars );
	}
	plot(
		av$log.m,
		av$log.vars,
		type="b",
		lty = 1,
		pch = "+",
		col = "royalblue",
		lwd = 2,
		main = main,
		xlab = xlab,
		ylab = ylab
	);
	#good <- av$m >= cutoffs[1] & av$m <= cutoffs[2]
	#points( av$log.m[good], av$log.vars[good], pch = "b", col = "royalblue", lwd = 2 ); # draw fitted points
	if ( plot.refline )
	{
		#abline(av$log.lm$coefficients[1]-av$log.lm$coefficients[1]*0.05, -1, col = "black", lty = 2 ); # reference line (5% of intercept below) (useful when both reference and fitting line are drawn)
		abline(av$log.lm$coefficients[1], -1, col = "black", lty = 2, lwd = 1 ); # reference line (5% of intercept below)
	}
	if ( plot.fitline )
	{
		abline(av$log.lm, col = "darkred", lty = 2, lwd = 1 ); # fitting line
	}
	if ( plot.cutlines )
	{
		abline(v=log10(av$cutoffs), col = "gray70", lty = 2, lwd = 1 ); # cut-off lines
	}
	grid();

	return( invisible(av) );
}

## MG_PLOT_LRD.RS
##
mg_plot_lrd.rs <- function( rs, n.levels = NULL, min.points = 3, cutoffs = 10^c(0.7, 3.5), main = NULL, xlab = "m", ylab = "R/S", col = "royalblue", plot.fitline = FALSE, plot.cutlines = TRUE )
{
	if (
		!is( rs, "mg_lrd_rs" )
	)
	{
		rs <- as.vector(rs);
		rs <- mg_lrd_rs(
			rs,
			n.levels = ifelse( !is.null(n.levels), n.levels, log10(length(rs)) ),
			min.points = min.points,
			cutoffs = cutoffs,
			ret.m = TRUE,
			ret.rs = TRUE,
			ret.log.m = TRUE,
			ret.log.rs = TRUE
		);
	}
	if (
		( is.null(rs$m) && is.null(rs$rs) )
		&& ( is.null(rs$log.m) && is.null(rs$log.rs) )
	)
	{
		stop( "[MG::PLOT::LRD::RS] Cannot plot R/S without (x,y) values." );
	}

	if ( is.null( rs$log.m ) )
	{
		rs$log.m <- log10( rs$m );
		rs$log.rs <- log10( rs$rs );
	}
	plot(
		rs$log.m,
		rs$log.rs,
		type="b",
		lty = 1,
		pch = "+",
		col = "royalblue",
		lwd = 2,
		main = main,
		xlab = xlab,
		ylab = ylab
	);
	#good <- rs$m >= cutoffs[1] & rs$m <= cutoffs[2]
	#points( rs$log.m[good], rs$log.rs[good], pch = "b", col = "royalblue", lwd = 2 ); # draw fitted points

#NOTE: To plot a refline we must know the real H (=d+1/2) value!!
#	if ( plot.refline )
#	{
#		#abline(rs$log.lm$coefficients[1]-rs$log.lm$coefficients[1]*0.05, -1, col = "black", lty = 2 ); # reference line (5% of intercept below) (useful when both reference and fitting line are drawn)
#		abline(rs$log.lm$coefficients[1], rs$H, col = "black", lty = 2, lwd = 1 ); # reference line (5% of intercept below)
#	}
	if ( plot.fitline )
	{
		abline(rs$log.lm, col = "darkred", lty = 2, lwd = 1 ); # fitting line
	}
	if ( plot.cutlines )
	{
		abline(v=log10(rs$cutoffs), col = "gray70", lty = 2, lwd = 1 ); # cut-off lines
	}
	grid();

	return( invisible(rs) );
}

## MG_PLOT_LRD.PERIODOGRAM
##
mg_plot_lrd.periodogram <- function( pgram, cutoff = 0.10, method = c("std", "cum"), main = NULL, xlab = "frequency", ylab = "periodogram", col = "royalblue", plot.fitline = FALSE, plot.cutlines = TRUE )
{
	if ( !is( pgram, "mg_lrd_periodogram" ) )
	{
		method <- match.arg( method );

		pgram <- as.vector(pgram);
		if ( method == "std" )
		{
			pgram <- mg_lrd_periodogram(
				pgram,
				cutoff = cutoff,
				ret.freqs = TRUE,
				ret.pgrams = TRUE,
				ret.log.freqs = TRUE,
				ret.log.pgrams = TRUE
			);
		} else {
			pgram <- mg_lrd_periodogram.cum(
				pgram,
				cutoff = cutoff,
				ret.freqs = TRUE,
				ret.pgrams = TRUE,
				ret.log.freqs = TRUE,
				ret.log.pgrams = TRUE
			);
		}
	}
	if (
		( is.null(pgram$freqs) && is.null(pgram$pgrams) )
		&& ( is.null(pgram$log.freqs) && is.null(pgram$log.pgrams) )
	)
	{
		stop( "[MG::PLOT::LRD::PERIODOGRAM] Cannot plot Periodograms without (x,y) values." );
	}

	if ( is.null( pgram$log.freqs ) )
	{
		pgram$log.freqs <- log10( pgram$freqs );
		pgram$log.pgrams <- log10( pgram$pgrams );
	}
	n <- length( pgram$log.freqs );
	# Plots used points
	plot(
		pgram$log.freqs[1:(n * pgram$cutoff)],
		pgram$log.pgrams[1:(n * pgram$cutoff)],
		type="p",
		lty = 1,
		pch = "+",
		col = "royalblue",
		lwd = 1,
		main = main,
		xlab = xlab,
		ylab = ylab
	);
	# Plots unused points
	points(
		pgram$log.freqs[(n * pgram$cutoff):n],
		pgram$log.pgrams[(n * pgram$cutoff):n],
		pch = "+",
		col = "gray70",
		lwd = 1
	);
	#good <- rs$m >= cutoffs[1] & rs$m <= cutoffs[2]
	#points( rs$log.m[good], rs$log.rs[good], pch = "b", col = "royalblue", lwd = 2 ); # draw fitted points

#NOTE: To plot a refline we must know the real H (=d+1/2) value!!
#	if ( plot.refline )
#	{
#		#abline(rs$log.lm$coefficients[1]-rs$log.lm$coefficients[1]*0.05, -1, col = "black", lty = 2 ); # reference line (5% of intercept below) (useful when both reference and fitting line are drawn)
#		abline(rs$log.lm$coefficients[1], rs$H, col = "black", lty = 2, lwd = 1 ); # reference line (5% of intercept below)
#	}
	if ( plot.fitline )
	{
		abline(pgram$log.lm, col = "darkred", lty = 3, lwd = 1 ); # fitting line
	}
	if ( plot.cutlines )
	{
		#abline(v=log10(pgram$freqs[n*pgram$cutoff]), col = "gray70", lty = 2, lwd = 1 ); # cut-off lines
		abline( v = pgram$log.freqs[n*pgram$cutoff], col = "gray70", lty = 2, lwd = 1 ); # cut-off lines
	}
	grid();

	return( invisible(pgram) );
}
