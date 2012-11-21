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

####
## Generic functions.
##
##  Marco Guazzone (marco.guazzone@gmail.com)
##
####

## Plots the empirical CDF of x.
## Params
##   * x: data values.
##   * title: plot title.
##   * xlabel: label of x-axe
##   * ylabel: label of y-axe
##   * col: foreground color.
mg_ecdfPlot <- function( x, title = "Empirical CDF", xlabel = "x", ylabel = "Fn(x)", col = "red" )
{
	Fn <- ecdf( x )
	plot(
		Fn,
		verticals = FALSE,
		#col.p = 'blue',
		do.points = FALSE,
		#digits = 5,
		col.h = col,
		col.v = col,
		lwd = 2,
		#col.v='bisque',
		main = title,
		xlab = xlabel,
		ylab = ylabel
	)
	grid()
}

## Plots the empirical CCDF of x.
## Params
##   * x: data values.
##   * title: plot title.
##   * xlabel: label of x-axe
##   * ylabel: label of y-axe
##   * col: foreground color.
mg_ccdfPlot <- function( x, title = "Empirical CCDF", xlabel = "x", ylabel = "Fn(x)}", col = "red" )
{
	x <-sort(x)
	Gn <- ecdf( x )
	plot(
		x,
		1-Gn(x),
		col = col,
		type = 'l',
		lwd = 2,
		main = title,
		xlab = xlabel,
		ylab = ylabel
	)
	grid()
}

## Plots the empirical LL-CCDF of x.
## Params
##   * x: data values.
##   * title: plot title.
##   * xlabel: label of x-axe
##   * ylabel: label of y-axe
##   * col: foreground color.
mg_llcdPlot <- function( x, title = "Empirical Log-Log CCDF", xlabel = "x", ylabel = "\bar{Fn(x)}", col = "red" )
{
	x <-sort(x)
	Fn <- ecdf( x )
	plot(
		x,
		1-Fn(x),
		log = 'xy',
		col = col,
		type = 'l',
		lwd = 2,
		main = title,
		xlab = xlabel,
		ylab = ylabel
	)
	grid()
}

## Plots a frequency histogram (i.e. a frequency bar plot)
mg_freqHistPlot <- function( x, breakVals = "Sturges", title = "Histogram", xlabel = "x", ylabel = "Frequency", color = "red" )
{
	mg_histPlot(
		x,
		breakVals = breakVals,
		showFreq = TRUE,
		col = color,
		main = title,
		xlab = xlabel,
		ylab = ylabel
	)
#	hist(
#		x,
#		breaks = breakVals,
#		freq = TRUE,
#		col = color,
#		main = title,
#		xlab = xlabel,
#		ylab = ylabel
#	)
#	grid()
}

## Plots a density histogram (i.e. a frequency bar plot)
mg_densHistPlot <- function( x, breakVals = "Sturges", title = "Histogram", xlabel = "x", ylabel = "Frequency", color = "red" )
{
	mg_histPlot(
		x,
		breakVals,
		showFreq,
		color,
		title,
		xlabel,
		ylabel
	);
#	hist(
#		x,
#		breaks = breakVals,
#		freq = FALSE,
#		col = color,
#		main = title,
#		xlab = xlabel,
#		ylab = ylabel
#	)
#	grid()
}

## Plots a histogram (i.e. a frequency bar plot)
## Params
##   * x: data values.
##   * title: plot title.
##   * xlabel: label of x-axe
##   * ylabel: label of y-axe
##   * col: foreground color.
mg_histPlot <- function( x, breakVals = "Sturges", showFreq = TRUE, title = "Histogram", xlabel = "x", ylabel = "Frequency", col = "red" )
{
	hist(
		x,
		breaks = breakVals,
		freq = showFreq,
		col = colors,
		main = title,
		xlab = xlabel,
		ylab = ylabel
	)
	grid()
}

## Plots a trend (i.e. a time series): data againts t)
## Params
##   * t: time values.
##   * data: data values.
##   * title: plot title.
##   * xlabel: label of x-axe
##   * ylabel: label of y-axe
##   * col: foreground color.
mg_trendPlot <- function( t, data, title = "Trend", xlabel = "t", ylabel = "x", col = "red" )
{
	jan1970 <- ISOdatetime( 1970, 1, 1, 0, 0, 0,'GMT' );
	plot(
		t + jan1970,
		data,
		#type = 'h',
		type = 'l',
		col = col,
		main = title,
		xlab = xlabel,
		ylab = ylabel,
		xaxt = 'n'
	);
	r <- as.POSIXct(range(t)+jan1970)
	axis.POSIXct( 1, at = seq( r[1], r[2], by = 'month' ), format = "%b" )
	grid()
}
