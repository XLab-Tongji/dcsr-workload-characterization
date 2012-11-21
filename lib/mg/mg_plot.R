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

## MG_PLOT
##
## SUMMARY
##  Generic plotting functions.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

source( "lib/mg/mg_eda.R" );
source( "lib/mg/mg_fit.R" );
source( "lib/mg/mg_gof_qq.R" );
source( "lib/mg/mg_gof_pp.R" );
source( "lib/mg/mg_plot_lrd.R" );

## MG_PLOT_FREQHIST(x,breaks,title,xlabel,ylabel,col,border.col,show.polygon,polygon.col,show.rug)
##
## SUMMARY
##  Plots a frequency histogram (i.e. a frequency bar plot)
##
mg_plot_freqHist <- function( x, breaks = "Sturges", title = "Histogram", xlabel = "x", ylabel = "Frequency", col = "gray90", border.col = "gray", show.polygon = TRUE, polygon.col = "royalblue", show.rug = TRUE )
{
#	h <- mg_plot_hist(
#		x,
#		breaks = breaks,
#		showFreq = TRUE,
#		col = col,
#		border.col = border.col,
#		title = title,
#		xlab = xlabel,
#		ylab = ylabel,
#		show.rug = show.rug
#	)
	if ( show.polygon )
	{
		# We must calculate the Y range for the polygonal

		h <- hist( x, breaks = breaks, plot = FALSE );

		# compute the frequency polygon
		diffBreaks <- h$mids[2] - h$mids[1]
		xx <- c( h$mids[1]-diffBreaks, h$mids, tail(h$mids,1)+diffBreaks )
		#yy <- c(0, h$density, 0)
		yy <- c(0, h$counts, 0)

		# draw the histogram 
		hist(
			x,
			breaks = breaks,
			freq = TRUE,
			xlim = range(xx),
			border = border.col,
			col = col,
			main = title,
			xlab = xlabel,
			ylab = ylabel
		);

		# adds the frequency polygon
		lines( xx, yy, lwd = 2, col = polygon.col );
	} else {
		h <- hist(
			x,
			breaks = breaks,
			freq = TRUE,
			#xlim = range(x),
			#border = border.col,
			col = col,
			main = title,
			xlab = xlabel,
			ylab = ylabel
		);
	}
	if ( show.rug )
	{
		rug( x );
	}
#	hist(
#		x,
#		breaks = breaks,
#		freq = TRUE,
#		col = col,
#		main = title,
#		xlab = xlabel,
#		ylab = ylabel
#	)
	grid();

	return( invisible( h ) );
}

## MG_PLOT_DENSITY(x,breaks,title,xlabel,ylabel,col,border.col,plot.density,density.kernel,density.col,show.rug)
##
## SUMMARY
##  Plots a density histogram (i.e. a frequency bar plot)
##
mg_plot_densHist <- function( x, breaks = "Sturges", title = "Histogram", xlabel = "x", ylabel = "Density", col = "gray90", border.col = "gray", plot.density = TRUE, density.kernel = "gaussian", density.width = NULL, density.col = "royalblue", show.rug = TRUE )
{
#	h <- mg_plot_hist(
#		x,
#		breaks = breaks,
#		showFreq = FALSE,
#		col = col,
#		border.col = border.col,
#		title = title,
#		xlab = xlabel,
#		ylab = ylabel,
#		show.rug = show.rug
#	);
	if ( plot.density )
	{
		# We must calculate the Y range for the histogram/density plot

		d <- density( x, kernel = density.kernel, width = density.width ); #FIXME: what default kernel shall I use?
		#lines( density( x ), lwd="2", col=density.col );

		# draw the histogram 
		h <- hist(
			x,
			prob = TRUE,
			breaks = breaks,
			#xlim = range( d$x ),
			ylim = range( d$y ),
			border = border.col,
			col = col,
			main = title,
			xlab = xlabel,
			ylab = ylabel
		);

		# adds the frequency polygon
		lines( d, lwd = 2, col = density.col );
	} else {
		# draw the histogram 
		h <- hist(
			x,
			prob = TRUE,
			breaks = breaks,
			border = border.col,
			col = col,
			main = title,
			xlab = xlabel,
			ylab = ylabel
		);
	}
	if ( show.rug )
	{
		rug( x );
	}
#	hist(
#		x,
#		breaks = breaks,
#		freq = FALSE,
#		col = col,
#		main = title,
#		xlab = xlabel,
#		ylab = ylabel
#	)
	grid();

	return( invisible( h ) );
}

## MG_PLOT_HIST(x,breaks,show.freq,title,xlabel,ylabel,col,border.col,show.rug)
##
## SUMMARY
##  Plots a histogram (i.e. a frequency bar plot)
##
## PARAMS
##  * x: data values.
##  * title: plot title.
##  * xlabel: label of x-axe
##  * ylabel: label of y-axe
##  * col: foreground color.
##
mg_plot_hist <- function( x, breaks = "Sturges", show.freq = TRUE, title = "Histogram", xlabel = "x", ylabel = "Frequency", col = "gray90", border.col = "gray", show.rug = TRUE )
{
	h <- hist(
		x,
		breaks = breaks,
		freq = show.freq,
		plot = TRUE,
		col = col,
		border = border.col,
		main = title,
		xlab = xlabel,
		ylab = ylabel
	)
	if ( show.rug )
	{
		rug( x );
	}

	return( invisible( h ) );
}

mg_plot_bar <- function( lengths, bar.names, show.polygon = TRUE, horiz = FALSE, ... )
{
	res <- barplot(
		lengths,
		names.arg = bar.names,
		col = "gray90",
		border = "gray",
		horiz = horiz,
		...
	);
	if ( show.polygon )
	{
		diff.breaks <- res[2] - res[1];
		if ( horiz )
		{
			lines(
				c( 0, lengths, 0 ),
				c( res[1] - diff.breaks, res, tail(res,1) + diff.breaks ),
				col = "royalblue",
				lwd = 2
			);
		} else {
			lines(
				c( res[1] - diff.breaks, res, tail(res,1) + diff.breaks ),
				c( 0, lengths, 0 ),
				col = "royalblue",
				lwd = 2
			);
		}
	}
	grid();
	return( invisible( res ) );
}

## MG_PLOT_TREND(t,data,title,xlabel,ylabel,col)
##
## SUMMARY
##  Plots a trend (i.e. a time series): data againts t)
##
## PARAMS
##  * t: time values.
##  * data: data values.
##  * title: plot title.
##  * xlabel: label of x-axe
##  * ylabel: label of y-axe
##  * col: foreground color.
##
mg_plot_trend <- function( data, title = "", xlabel = "index", ylabel = "observation", col = "royalblue" )
{
	plot(
		data,
		type = 'h',
		#type = 'l',
		col = col,
		main = title,
		xlab = xlabel,
		ylab = ylabel
	);

	grid();
}

## MG_PLOT_TREND(t,data,title,xlabel,ylabel,col)
##
## SUMMARY
##  Plots a trend (i.e. a time series): data againts t)
##
## PARAMS
##  * t: time values.
##  * data: data values.
##  * title: plot title.
##  * xlabel: label of x-axe
##  * ylabel: label of y-axe
##  * col: foreground color.
##
mg_plot_trend.ts <- function( t, data, title = "Trend", xlabel = "time", ylabel = "x", col = "royalblue" )
{
	jan1970 <- ISOdatetime( 1970, 1, 1, 0, 0, 0,"GMT" );
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
	r <- as.POSIXct(range(t)+jan1970);
	axis.POSIXct( 1, at = seq( r[1], r[2], by = "month" ), format = "%b" );
	grid();
}

## MG_PLOT_LAG(x,lag=1,title,xlab,ylab,col)
##
## SUMMARY
##  Draws a Lag-plot.
##
mg_plot_lag <- function( x, lag = 1, do.refline = TRUE, title = NULL, xlab = paste("lag",lag), ylabel = "x", col = "royalblue", col.refline = "gray70" )
{
	lag.plot(
		x,
		lags = lag,
		diag = do.refline,
		diag.col = col.refline,
		col = col,
		pch = "+"
	);
	grid();
}

## MG_PLOT_HILL(x,conf.level,title,xlabel,ylabel,col)
##
## SUMMARY
##  Plots a Hill plot.
##
## PARAMS
##  * x: data values.
##  * conf.level: cofidence level.
##  * title: plot title.
##  * xlabel: x-axis label.
##  * ylabel: y-axis label.
##  * col: line color.
##
mg_plot_hill <- function( x, conf.level = 0.95, title = "", xlabel = "Tail size", ylabel = "Hill estimate", col = "royalblue" )
{
	n <- length( x );
#	y <- c();
##n<-1000
#	for ( k in 1:(n-1) )
#	{
#		e <- mg_fit_hillEstimator( x, k );
#		y[k] <- e["shape"];
#	}
	e <- mg_fit_hillEstimator( x, n-1, conf.level );
	k <- 1:e$k.max;
	ylim <- NULL;
	env.draw <- FALSE; # draw envelope?
	if ( !all( is.na( e$conf.int$lower ) ) && !all( is.na( e$conf.int$upper ) ) )
	{
		ylim <- c( min( e$conf.int$lower ), max( e$conf.int$upper ) );
		env.draw <- TRUE;
	}
	plot(
		#1:(n-1),
		#y,
		k,
		e$statistic,
		type = "l",
		lty = "solid",
		lwd = 2,
		col = col,
		main = title,
		xlab = xlabel,
		ylab = ylabel,
		ylim = ylim
	);

	# Draws confidence envelope
#	if ( !all( is.na( e$conf.int$upper ) ) )
#	{
#		lines( k, e$conf.int$upper, lty = "dashed", col = "black", lwd = 1 );
#	}
#	if ( !all( is.na( e$conf.int$lower ) ) )
#	{
#		lines( k, e$conf.int$lower, lty = "dashed", col = "black", lwd = 1 );
#	}
	if ( env.draw )
	{
		lines( k, e$conf.int$upper, lty = "dashed", col = "darkred", lwd = 1 );
		lines( k, e$conf.int$lower, lty = "dashed", col = "darkred", lwd = 1 );
	}

	grid();

	return( invisible( e ) );
}

## MG_PLOT_PP
##
## SUMMARY
##  Plot a P-P plot.
##
## DESCRIPTION
##  The P-P (Probability-Probability) plot plots sorted F(x(i)) against (i-1/2)/n.
mg_plot_pp <- function( y, dist, ..., use.points = TRUE, show.ppline = TRUE, show.pppoints = TRUE, log  = "", col = "royalblue", xlab = paste(dist,"CDF",sep=" "), ylab = paste("Empirical CDF") )
{
	y <- sort(y);
	xx <- mg_dists_cdf( dist, y, ... );
	yy <- ppoints(xx)[order(order(y))]; # ((1:n)-0.5)/n

	#plot( 0,0, type="n" );
	plot(
		xx,
		yy,
		col = col,
		lwd = ifelse( use.points, 1, 2),
		pch = "+",
		type = ifelse( use.points, "p", "l" ),
		xlab = xlab,
		ylab = ylab,
		log = log
	);

        #Axis( x=c(0,1), side=1 );
        #Axis( x=c(0,1), side=2 );

	#ppstat <- list();
	#p1 <- .25;
	#p3 <- .75;
	#data.qq <- quantile( xx, probs = c( p1, p3 ) );
	#distr.qq <- mg_dists_invcdf( dist, xx, ... );
	#ppstat$q1 = c( distr.qq[1], data.qq[1] );
	#ppstat$q3 = c( distr.qq[2], data.qq[2] );
	## Draw a line passing through 25% and 75% percentile.
	##   slope = (y2-y1)/(x2-x1)
	##   int   = y1 - x1*(y2-y1)/(x2-x1)
	#ppstat$slope <- (data.qq[2] - data.qq[1]) / (distr.qq[2]-distr.qq[1]);
	#ppstat$int <- data.qq[1] - distr.qq[1]*ppstat$slope;
	#linefun <- function(x) { return( ppstat$int + ppstat$slope * x ); };
	ppstats <- mg_gof_pp.stat2( xx, yy );

	if ( show.ppline )
	{
		abline(ppstats$int, ppstats$slope, col = "black", lwd = 1, lty = "dashed" );
		##abline( .25, ppstat$slope, col = "magenta", lwd = 1, lty = "dashed" );
		##l <- lm( yy ~ xx );
		##abline( l$coefficients, col = "magenta", lwd = 1, lty = "dashed" );
		#abline( 0, 1, col = "red", lwd = 1, lty = "dashed" );
		##abline( qqstat$int, 1, col = "magenta", lwd = 1, lty = "dotted", ... );
	}
	if ( show.pppoints )
	{
		#points( c(ppstats$p1[1], ppstats$p3[1]), c(ppstats$p1[2], ppstats$p3[2]), pch = 20, cex = 2, col = "red", ... );
		points( c(ppstats$p1[1], ppstats$p3[1]), c(ppstats$p1[2], ppstats$p3[2]), pch = 20, col = "red" );
		#points( c(linefun(ppstat$q1[1]), linefun(ppstat$q3[1])), c(linefun(ppstat$q1[2]), linefun(ppstat$q3[2])), pch = 20, col = "red" );
	}

	grid();

	return( invisible( ppstats ) );
}

## MG_PLOT_QQ(x,dist,...,use.points,show.qqline,shoq.qqpoints,col)
##
## SUMMARY
##  Plot a One-sample QQ-plot.
##
mg_plot_qq <- function( x, dist, ..., show.qqline = TRUE, show.qqpoints = TRUE, col = "royalblue", log = "" )
{
	q <- try( mg_dists_invcdf( dist, ppoints(x), ... ) );
	if ( is( q, "try-error" ) )
	{
		stop( "[MG::PLOT::QQ] Error while calculating theoretical quantiles for distribution '", dist, "': ", geterrmessage() );
	}

	return(
		mg_plot_qq2(
			q[order(order(x))], # order(oder(x)) == rank(x)
			#mg_dists_invcdf( dist, ppoints(x), ... ),
			x,
			show.qqline = show.qqline,
			show.qqpoints = show.qqpoints,
			col = col,
			xlab = dist,
			ylab = "x",
			log = log
		)
	);
}

## MG_PLOT_QQ2(x,y,use.points,show.qqline,shoq.qqpoints,col,...)
##
## SUMMARY
##  Plot a Two-samples QQ-plot.
##
mg_plot_qq2 <- function( x, y, use.points = TRUE, show.qqline = TRUE, show.qqpoints = TRUE, log = "", col = "royalblue", xlab = "x", ylab = "y", ... )
{
	# Sanitizes data
	good <- !is.na(x); # index set of NA values for x
	if ( length(good) < length(x) )
	{
		# Removes NA entries
		x <- x[good];
		y <- y[good];
	}
	good <- !is.na(y) # idex set of NA values for y
	if ( length(good) < length(y) )
	{
		# Removes NA entries
		x <- x[good];
		y <- y[good];
	}

	# sanity check
	if ( length(x) == 0 || length(y) == 0 )
	{
		return( NULL );
	}

	if (length(log) && log != "") {
		log <- strsplit(log, NULL)[[1]]
		if ("x" %in% log )
		{
			if ( any(ii <- x <= 0 & !is.na(x)))
			{
				n <- as.integer(sum(ii))
				warning(sprintf(ngettext(n, "%d x value <= 0 omitted from logarithmic plot", "%d x values <= 0 omitted from logarithmic plot"), n), domain = NA)
				x[ii] <- NA
			}
			x <- log(x);
		}
		if ("y" %in% log )
		{
			if ( any(ii <- y <= 0 & !is.na(y)))
			{
				n <- as.integer(sum(ii))
				warning(sprintf(ngettext(n, "%d y value <= 0 omitted from logarithmic plot", "%d y values <= 0 omitted from logarithmic plot"), n), domain = NA)
				y[ii] <- NA
			}
			y <- log(y);
		}
	}

	qqp <- qqplot(
		x,
		y,
		col = col,
		lwd = ifelse( use.points, 1, 2),
		pch = "+",
		type = ifelse( use.points, "p", "l" ),
		xlab = xlab,
		ylab = ylab,
		...
	);

#	qqslope <- mg_plot_qqline(
#		qqp$x,
#		qqp$y,
#		plot.it = show.qqline,
#		show.qqpoints = show.qqpoints
#	);
#	#qqslope <- qqslope$slope;
##	if ( show.qqline )
##	{
##		mg_plot_qqline( x, y, show.qqpoints = show.qqpoints );
##	}
##	else if ( show.qqpoints )
#	if ( !show.qqline && show.qqpoints )
#	{
#		mg_plot_qqpoints( qqp$x, qqp$y );
#	}

	qqstat <- mg_gof_qq.stat2( qqp$x, qqp$y );

	if ( show.qqline )
	{
		abline(qqstat$int, qqstat$slope, col = "black", lwd = 1, lty = "dashed", ... );
		#abline( qqstat$int, 1, col = "magenta", lwd = 1, lty = "dotted", ... );
	}
	if ( show.qqpoints )
	{
		#points( c(qqstat$q1[1], qqstat$q3[1]), c(qqstat$q1[2], qqstat$q3[2]), pch = 20, cex = 2, col = "red", ... );
		points( c(qqstat$q1[1], qqstat$q3[1]), c(qqstat$q1[2], qqstat$q3[2]), pch = 20, col = "red", ... );
	}

	grid();

	return( invisible( qqstat ) );
}

# confidence.band (Derek Young and David Hunter)
# based on Ellipses, by J. Fox and G. Monette, from
# car package

#confidence.band = function(model, levels=0.95, segments=50, col.points=palette()[1], 
#                           col.line=palette()[1], col.bands=palette()[2], 
#                           lty.line=1, lty.bands=2, ...) {
#  if (attr(model$terms,"intercept")!=1 || length(model$coef) !=2) {
#    stop(paste("condifence.bands only works for simple linear regression\n",
#               "with one predictor and an intercept"))
#  }
#  plot(model$model[,2:1], col=col.points, ...)
#  abline(model, col=col.line, lty=lty.line, lwd=2)
#  angles=(0:segments)*pi/segments
#  halfcircle = cbind(cos(angles), sin(angles))  
#  chol.shape = chol(vcov(model))
#  slopes = (halfcircle %*% chol.shape)[,2]
#  angles = angles+angles[which.max(slopes)]
#  halfcircle = cbind(cos(angles), sin(angles))  
#  center = model$coef
#  radius = sqrt(2*qf(levels, 2, df.residual(model)))
#  for (r in radius) {
#    for (i in 1:2) {
#      halfcircle = -halfcircle
#      ellipse = sweep(r*(halfcircle %*% chol.shape), 2, center, "+")
#      int = ellipse[,1]
#      slope = ellipse[,2]
#      x = -diff(int)/diff(slope)
#      y = int[-1]+slope[-1]*x
#      lines(x, y, lwd=2, lty=lty.bands, col=col.bands)
#    }
#  }
#}

#mg_plot_qqline <- function( x, y, col = "black", plot.it = TRUE, show.qqpoints = TRUE, qqpoints.col = "red", conf.level = 0.95, show.confint = TRUE, confint.col = "darkgreen", ... )
#{
#	# the quantile line
#	qx <- quantile(x[!is.na(x)], c(0.25, 0.75)); # Q1,Q3 for x
#	qy <- quantile(y[!is.na(y)], c(0.25, 0.75)); # Q1,Q3 for y
#	slope <- diff(qy)/diff(qx); # dy/dx
#	int <- qy[1] - slope * qx[1]; # intercept
#
#	# TODO: add a new type of line: robust regression line
#	# (note: x and y refers to arguments version)
#	# TO BE TESTED
#        #if ( !require("MASS") )
#	#{
#	#	stop("MASS package not available");
#	#}
#        #coef <-MASS::coefficients( MASS::rlm( y ~ x ) );
#        #int <-coef[1];
#        #slope <-coef[2];
#	#abline( int, slope );
#
#	# TODO: add a new type of line: 45-degree reference line
#	# TO BE TESTED
#	#abline( 0, 1 ):
#
#	# TODO: add confidence interval lines (the envelope)
#	# (note: x and y refers to arguments version)
#	# TO BE TESTED
##	se <- slope / sqrt(length(x));
##	qqz <- qnorm(1 - (1 - conf.level)/2);
##	qqv <- int + slope * x;
##	u <- qqv + se * qqz;
##	l <- qqv - se * qqz;
#
#	qqz <- qt( 1 - (1 - conf.level)/2, length(x) - 1 );
#	qqv <- int + slope * x;
#	s <- sd( qqv );
#	m <- mean( qqv );
#	se <- sqrt( s ) / sqrt( length(qqv) );
#	l <- qqv + se * qqz;
#	u <- qqv - se * qqz;
#
##	l <- x - se
##	P <- ppoints(length(x))
##	zz<-qnorm(1-(1-conf.level)/2)
##	se<-(slope/d.function(x, ...))*sqrt(P*(1-P)/n)
##	fit.value<-a+b*z
##	upper<-fit.value+zz*Se
##	lower<-fit.value-zz*SE
##	lines(z, upper, lty=2, lwd=lwd/2, col=col)
##	lines(z, lower, lty=2, lwd=lwd/2, col=col)
#
##	centerx <- (qx[1] + qx[2])/2;
##	centery <- (qy[1] + qy[2])/2;
##	maxx <- max(x);
##	minx <- min(x);
##	maxy <- centery + slope*(maxx - centerx);
##	miny <- centery - slope*(centerx - minx);
##
##	mx <- c(minx, maxx);
##	my <- c(miny, maxy);
#
#	if ( plot.it )
#	{
#		abline(int, slope, col = col, lwd = 1, lty = "dashed", ... );
#
#		if ( show.qqpoints )
#		{
#			points( qx, qy, pch = 20, cex = 2, col = qqpoints.col );
#		}
#		if ( show.confint )
#		{
#			lines( x, u, lty = "dotted", lwd = 1, col = confint.col );
#			lines( x, l, lty = "dotted", lwd = 1, col = confint.col );
#		}
##		lines( qx, qy, pch = 20, cex = 2, col = "yellow" );
#
#		grid();
#	}
#
#	return( invisible( list( slope = slope, qqx = x, qqy = y ) ) );
#}
#
#mg_plot_qqpoints <- function( x, y, col = "red", pch = 20, ... )
#{
#	x <- quantile(x[!is.na(x)], c(0.25, 0.75));
#	y <- quantile(y[!is.na(y)], c(0.25, 0.75));
#
#	points( x, y, pch = pch, cex = 2, col = col );
#
#	grid();
#
#	return( invisible( list( qqx = x, qqy = y ) ) );
#}

## MG_PLOT_ECDF(x,...,ylab,verticals,do.points,log,col.line,col.points,col.01line)
##
## SUMMARY
##  Plots the empirical CDF of x.
##
## PARAMS
##  * x: data values.
##
mg_plot_ecdf <- function ( x, ..., ylab = "Fn(x)", verticals = FALSE, do.points = TRUE, log = "", col.line = "royalblue", col.points = "royalblue", col.01line = "gray70", add = FALSE ) 
{
	if ( !is.stepfun(x) )
	{
		if ( is.numeric(x) )
		{
			sarg <- substitute(x);
			x <- mg_eda_ecdf(x);
			attr(x, "call") <- call("mg_eda_ecdf", sarg);
		} else {
			stop("[MG::PLOT::ECDF] called with wrong type of argument 'x'");
		}
	}
	if ( !is.null(log) && log != "" )
	{
		knF <- knots( x );
		knF <- knF[knF > 0];
		Fn <- x(knF);

		if ( add )
		{
			# Note: don't pass log option since assume the axis is already in log scale

			# (taken from plot.stepfun)
			rx <- range(knF)
			dr <- max(0.08 * diff(rx), median(diff(knF)))
			xlim <- rx + dr * c(-1, 1)
			knF <- knF[xlim[1] - dr <= knF & knF <= xlim[2] + dr]
			ti <- c(xlim[1] - dr, knF, xlim[2] + dr)
			ti.l <- ti[-length(ti)]
			ti.r <- ti[-1]
			y <- x(0.5 * (ti.l + ti.r))
			n <- length(y)
			#Fn.kn <- x(knF)
			segments(ti.l, y, ti.r, y, col = col.line,  ylab = ylab,, ...)
		} else {
			plot( knF, Fn, type = "s", col = col.line, log = log, ylab = ylab, ... );
		}

		if ( verticals )
		{
			#TODO: how do we do it?
		}
		if ( do.points )
		{
			points( knF, Fn, col = col.points );
		}
		#plot.stepfun( knF, Fn, col.hor = col.line, col.vert = col.line, col.points = col.points, do.points = do.points, verticals = verticals, ylab = ylab, add = add, ... );
	} else {
		plot.stepfun(
			x,
			#...,
			ylab = ylab,
			verticals = verticals,
			do.points = do.points,
			col.points = col.points,
			col.hor = col.line,
			col.vert = col.line,
			add = add,
			...
		);
		abline( h = c(0, 1), col = col.01line, lty = 2 );
	}
	grid();
}

## MG_PLOT_ECCDF(x,...,ylab,verticals,do.points,log,col.line,col.points,col.01line)
##
## SUMMARY
##  Plots the empirical CCDF of x.
##
## PARAMS
##  * x: data values.
##
mg_plot_eccdf <- function ( x, ..., ylab = "1-Fn(x)", verticals = FALSE, do.points = TRUE, log = "", col.points = "royalblue", col.line = "royalblue", col.01line = "gray70", add = FALSE ) 
{
	if ( !is.stepfun(x) )
	{
		if ( is.numeric(x) )
		{
			sarg <- substitute(x);
			x <- mg_eda_eccdf(x);
			attr(x, "call") <- call("mg_eda_eccdf", sarg);
		} else {
			stop("[MG::PLOT::ECCDF] called with wrong type of argument 'x'");
		}
	}
	if ( !is.null(log) && log != "" )
	{
		knF <- knots( x );
		knF <- knF[knF > 0];
		Fn <- x(knF);

		if ( add )
		{
			# Note: don't pass log option since assume the axis is already in log scale

			# (taken from plot.stepfun)
			rx <- range(knF)
			dr <- max(0.08 * diff(rx), median(diff(knF)))
			xlim <- rx + dr * c(-1, 1)
			knF <- knF[xlim[1] - dr <= knF & knF <= xlim[2] + dr]
			ti <- c(xlim[1] - dr, knF, xlim[2] + dr)
			ti.l <- ti[-length(ti)]
			ti.r <- ti[-1]
			y <- x(0.5 * (ti.l + ti.r))
			n <- length(y)
			#Fn.kn <- x(knF)
			segments(ti.l, y, ti.r, y, col = col.line,  ylab = ylab,, ...)
		} else {
			plot( knF, Fn, type = "s", col = col.line, log = log, ylab = ylab, ... );
		}

		if ( verticals )
		{
			#TODO: how do we do it?
		}
		if ( do.points )
		{
			points( knF, Fn, col = col.points );
		}
	} else {
		plot.stepfun(
			x,
			...,
			ylab = ylab,
			verticals = verticals,
			do.points = do.points,
			col.points = col.points,
			col.hor = col.line,
			col.vert = col.line,
			add = add
		);
		abline( h = c(0, 1), col = col.01line, lty = 2 );
	}
	grid();
}

mg_plot_boxplot <- function( x, ..., col.fill = "gray90", col.border = "royalblue", title = NULL, names = NULL )
{
	boxplot(x,...,col=col.fill,border=col.border,pch="+", main=title, names = names);
	grid();

	return( invisible() );
}

mg_plot_aggregation.trend <- function( x, scale.start = 1, scale.end = floor(log10(length(x))), scale.incr = 0, scale.mult = 10, fit.method = c("none", "accumulate", "enlarge"), main = "", xlab = "index", ylab = "observation", plot.compact = FALSE )
{
	fit.method <- match.arg(fit.method);
	scale.start <- max( as.integer(scale.start), 1 );
	scale.end <- max( as.integer(scale.end), scale.start );
	scale.incr <- max( as.integer(scale.incr), 0 );
	scale.mult <- max( as.integer(scale.mult), 1 );

	scale.n <- ceiling( (scale.end-scale.start+1)/scale.mult );

	colors <- NULL;
	if ( plot.compact )
	{
		colors <- rainbow( scale.n );
	}
	else
	{
		colors <- rep( "royalblue", scale.n );
	}

	ltype <- 1;
	labels <- c();
	ltypes <- c();

	x.len <- length(x);

#	if ( scale.start == 1 )
#	{
#		plot( x, type = "h", col = colors[1], main = main, xlab = xlab, ylab = ylab );
#		labels <- c( paste( "x", scale.start, sep = "" ) );
#		ltypes <- c( ltype );
#	}

	oldpar <- NULL; # old graphical parameters
	if ( !plot.compact && scale.start != scale.end )
	{
		max.rows <- 6;
		cols <- ceiling( scale.n / max.rows );
		rows <- ceiling( scale.n / cols );
		oldpar <- par( mfrow = c( rows, cols ) );
	}

	scale <- scale.start;
	for ( i in 1:scale.n )
	{
		agg.len <- x.len / scale;
		agg <- mg_eda_aggregate( x, n.groups = floor( agg.len ), fit.method = fit.method ); 
		if ( i == 1 || !plot.compact )
		{
			plot( agg, type = "h", col = colors[1], lty = 1, main = main, xlab = xlab, ylab = ylab );

			if ( !plot.compact )
			{
				# Adds a grid and a title for each subplot
				title( main = paste( "x", scale, sep = "" ) );
				grid();
			}
		} else {
			lines( agg, type = "h", col = colors[i], lty = ltype );
		}

		agg <- NULL; # clean-up mem

		labels <- append( labels, paste( "x", scale, sep = "" ) );
		ltypes <- append( ltypes, ltype );

		#scale <- scale.start + i * scale.step;
		scale <- (i - 1 + scale.start) * scale.mult + scale.incr;

		ltype <- ltype + 1;
		if ( ( ltype %% 6 ) == 1 )
		{
			# avoid solid line (used for data)
			ltype <- ltype + 1;
		}
	}

	if ( plot.compact )
	{
		# Adds a global legend and grid

		legend( x = "topright", legend = labels, col = colors, lty = 1 );
		grid();
	}

	if ( !is.null(oldpar) )
	{
		par( oldpar );
	}

	return( invisible( agg ) );
}

mg_plot_aggregation.ecdf <- function( x, scale.start = 1, scale.end = floor(log10(length(x))), scale.incr = 0, scale.mult = 10, fit.method = c("none", "accumulate", "enlarge"), log = "", main = "", xlab = "index", ylab = "observation" )
{
	fit.method <- match.arg(fit.method);
	scale.start <- max( as.integer(scale.start), 1 );
	scale.end <- max( as.integer(scale.end), scale.start );
	#scale.step <- max( as.integer(scale.step), 1 );
	scale.incr <- max( as.integer(scale.incr), 0 );
	scale.mult <- max( as.integer(scale.mult), 1 );

	#scale.n <- ceiling( (scale.end-scale.start+1)/scale.step );
	scale.n <- ceiling( (scale.end-scale.start+1)/scale.mult );

	cols <- rainbow( scale.n );

	if ( log != "" )
	{
		x <- x[x > 0];
	}

	x.len <- length(x);

	ltype <- 1;
	ltypes <- c();
	labels <- c();

	# Plots scaled data
	scale <- scale.start;
	for ( i in 1:scale.n )
	{
		agg.len <- x.len / scale;
		agg <- mg_eda_aggregate( x, n.groups = floor( agg.len ), fit.method = fit.method );
		if ( i == 1 )
		{
#			mg_plot_ecdf( agg, add = FALSE, verticals = TRUE, do.points = FALSE, col.line = cols[i], col.points = cols[i], lty = ltype, log = log, main = main, xlab );
			mg_plot_ecdf( agg, verticals = TRUE, do.points = FALSE, col.line = cols[1], col.points = cols[1], lty = ltype, main = main, xlab = xlab, ylab = ylab, log = log );
		} else {
			mg_plot_ecdf( agg, add = TRUE, verticals = TRUE, do.points = FALSE, col.line = cols[i], col.points = cols[i], lty = ltype, log = log );
		}

		agg <- NULL; # clean-up mem

		labels <- append( labels, paste( "x", scale, sep = "" ) );
		ltypes <- append( ltypes, ltype );

		#scale <- scale.start + i * scale.step;
		scale <- (i - 1 + scale.start) * scale.mult + scale.incr;

		ltype <- ltype + 1;
		if ( ( ltype %% 6 ) == 1 )
		{
			# avoid solid line (used for data)
			ltype <- ltype + 1;
		}
	}

	legend( x = "bottomright", legend = labels, col = cols, lty = ltypes );

	grid();
}

mg_plot_aggregation.eccdf <- function( x, scale.start = 1, scale.end = floor(log10(length(x))), scale.incr = 0, scale.mult = 10, fit.method = c("none", "accumulate", "enlarge"),log = "", main = "", xlab = "index", ylab = "observation" )
{
	fit.method <- match.arg(fit.method);
	scale.start <- max( as.integer(scale.start), 1 );
	scale.end <- max( as.integer(scale.end), scale.start );
	#scale.step <- max( as.integer(scale.step), 1 );
	scale.incr <- max( as.integer(scale.incr), 0 );
	scale.mult <- max( as.integer(scale.mult), 1 );

	#scale.n <- ceiling( (scale.end-scale.start+1)/scale.step );
	scale.n <- ceiling( (scale.end-scale.start+1)/scale.mult );

	cols <- rainbow( scale.n );

	if ( log != "" )
	{
		x <- x[x > 0];
	}

	x.len <- length(x);

	ltype <- 1;
	labels <- c();
	ltypes <- c();

#	# Plots unscaled data
#	if ( scale.start == 1 )
#	{
#		mg_plot_eccdf( x, verticals = TRUE, do.points = FALSE, col.line = cols[1], col.points = cols[1], lty = ltype, main = main, xlab = xlab, ylab = ylab, log = log );
#		#labels <- c( "x1" );
#		labels <- c( paste( "x", scale.start, sep = "" ) );
#		ltypes <- c( ltype );
#	}

	# Plots scaled data
	scale <- scale.start;
	for ( i in 1:scale.n )
	{
		agg.len <- x.len / scale;
		agg <- mg_eda_aggregate( x, n.groups = floor( agg.len ), fit.method = fit.method ); 

		if ( i == 1 )
		{
			mg_plot_eccdf( agg, verticals = TRUE, do.points = FALSE, col.line = cols[1], col.points = cols[1], lty = ltype, main = main, xlab = xlab, ylab = ylab, log = log );
		} else {
			mg_plot_eccdf( agg, add = TRUE, verticals = TRUE, do.points = FALSE, col.line = cols[i], col.points = cols[i], lty = ltype, log = log );
		}

		agg <- NULL; # clean-up mem

		labels <- append( labels, paste( "x", scale, sep = "" ) );
		ltypes <- append( ltypes, ltype );

		#scale <- scale.start + i * scale.step;
		scale <- (i - 1 + scale.start) * scale.mult + scale.incr;

		ltype <- ltype + 1;
		if ( ( ltype %% 6 ) == 1 )
		{
			# avoid solid line (used for data)
			ltype <- ltype + 1;
		}
	}

	legend( x = "bottomleft", legend = labels, col = cols, lty = ltypes );

	grid();
}

mg_plot.lorenz <- function(x, title = NULL, xlab = "% count", ylab = "% mass", col.count = "royalblue", col.mass = "forestgreen" )
{
	l <- mg_eda.lorenz(x);
	plot( l$count*100, l$mass*100, type = "l", col = "royalblue", lwd = 2, main = title, xlab = xlab, ylab = ylab );
	abline( 0, 1, lty = 2, col = "black", lwd = 1 ); # the "perfect equality" line
	abline( h=0, v = 100, lty = 2, col = "black", lwd = 1 ); # the "perfect inequality" line
	grid();

	return( invisible( l ) );
}

mg_plot.masscount <- function(x, title = NULL, xlab = "count CDF", ylab = "mass CDF", col.count = "royalblue", col.mass = "forestgreen" )
{
	l <- mg_eda.lorenz(x);
	plot( l$x, l$count, type = "l", col = col.count, lwd = 2, main = title, xlab = xlab, ylab = ylab );
	lines( l$x, l$mass, col = col.mass, lwd = 2 );

	legend( "bottomright", legend = c( "count", "mass" ), col = c( col.count, col.mass ), lty = 1 );

	grid();

	return( invisible( l ) );
}

## MG_PLOT_MULTIPLOT( x, distrs, plot.type, title=NULL, xlabel=NULL, ylabel=NULL, colors=NULL )
##
## SUMMARY
##  Draw multiple plots on a same graphic window.
##
mg_plot_multiplot <- function( x, distrs, plot.type, title = NULL, xlabel = NULL, ylabel = NULL, colors = NULL )
{
	x <- na.omit(x); # remove NAs

	n.distrs <- length( distrs );

	# sanity check
	if ( length(x) == 0 || n.distrs == 0 )
	{
		warning( "[MG::PLOT::MULTIPLOT] Not enough data for plotting (", length(x), ",", n.distrs, ")" );
		return(invisible(NULL));
	}

	# adjust colors
	if ( is.null( colors ) || length( colors ) <= n.distrs )
	{
		colors <- rainbow( 1 + n.distrs );
	}

	# Line types vector used for legend
	ltypes <- c(1); # initialize to "solid" line type
	col.idx <- 1; # color index
	plot.distrs <- c(); # names of plot distributions
	leg.place <- NULL; # legend placement

	if ( plot.type == MG_CONSTS_CCDF_PLOT || plot.type == MG_CONSTS_LLCD_PLOT )
	{
		x <- sort( x );
		x <- x[ x > 0 ]; # In log-log scale values must be > 0

		# sanity check
		if ( length(x) == 0 )
		{
			warning( "[MG::PLOT::MULTIPLOT] Not enough data to plot '", plot.type, "'" );
			return(invisible(NULL));
		}

		# Plot ECCDF
		mg_plot_eccdf(
			x,
			log = ifelse( plot.type == MG_CONSTS_LLCD_PLOT, "xy", "" ),
			col.line = colors[ col.idx ],
			#type = 'l',
			lty = "solid", # => lty = 1
			lwd = 2,
			main = ifelse( is.null(title), ifelse( plot.type == MG_CONSTS_LLCD_PLOT, "log-log CCDF", "CCDF" ), title ),
			xlab = ifelse( is.null(xlabel), ifelse( plot.type == MG_CONSTS_LLCD_PLOT, "log(x)", "x" ), xlabel ),
			ylab = ifelse( is.null(ylabel), ifelse( plot.type == MG_CONSTS_LLCD_PLOT, "log(CCDF(x))", "CCDF(x)" ), ylabel ),
			verticals = TRUE,
			do.points = FALSE
		);

		ltype <- 1; # line type progess
		for ( distr in names( distrs ) )
		{
			# sanity check
			if ( is.null( distrs[[ distr ]] ) || length( distrs[[ distr ]] ) == 0 )
			{
				next;
			}

			# Updates line type progress
			ltype <- ltype + 1;
			if ( ( ltype %% 6 ) == 1 )
			{
				# avoid solid line (used for data)
				ltype <- ltype + 1;
			}

			# Updates color
			col.idx <- col.idx + 1; 

			# plot CDF
			res <- try(
				lines(
					x,
					1 - do.call( "mg_dists_cdf",  c( list(distr), list(x), distrs[[ distr ]] ) ),
					col = colors[ col.idx ],
					lty = ltype,
					lwd = 2
				)
			);
			if ( is( res, "try-error" ) )
			{
				stop( "[MG::PLOT::MULTIPLOT] Error while drawing plot '", plot.type, "': ", geterrmessage() );
			}

			# Updates legend infos
			ltypes <- append( ltypes, ltype );
			plot.distrs <- append( plot.distrs, distr );
		}
		leg.place <- "topright";
	} else if ( plot.type == MG_CONSTS_CDF_PLOT || plot.type == MG_CONSTS_LLCDF_PLOT ) {
		x <- sort( x );

		# sanity check
		if ( length(x) == 0 )
		{
			warning( "[MG::PLOT::MULTIPLOT] Not enough data to plot '", plot.type, "'" );
			return(invisible(NULL));
		}

		# Plot ECDF
		mg_plot_ecdf(
			mg_eda_ecdf(x),
			verticals = TRUE,
			do.points = FALSE,
			log = ifelse( plot.type == MG_CONSTS_LLCDF_PLOT, "xy", "" ),
			col.line = colors[ col.idx ],
			#type = 'l',
			lty = "solid", # => lty = 1
			lwd = 2,
			main = ifelse( is.null(title), ifelse( plot.type == MG_CONSTS_LLCDF_PLOT, "log-log CDF", "CDF" ), title ),
			xlab = ifelse( is.null(xlabel), ifelse( plot.type == MG_CONSTS_LLCDF_PLOT, "log(x)", "x" ), xlabel ),
			ylab = ifelse( is.null(ylabel), ifelse( plot.type == MG_CONSTS_LLCDF_PLOT, "log(CDF(x))", "CDF(x)" ), ylabel )
		);

		ltype <- 1; # line type progess
		for ( distr in names( distrs ) )
		{
			# sanity check
			if ( is.null( distrs[[ distr ]] ) || length( distrs[[ distr ]] ) == 0 )
			{
				next;
			}

			# Updates line type progress
			ltype <- ltype + 1;
			if ( ( ltype %% 6 ) == 1 )
			{
				# avoid solid line (used for data)
				ltype <- ltype + 1;
			}

			# Updates color
			col.idx <- col.idx + 1; 

			# plot CDF
			res <- try(
				lines(
					x,
					do.call( "mg_dists_cdf",  c( list(distr), list(x), distrs[[ distr ]] ) ),
					col = colors[ col.idx ],
					lty = ltype,
					lwd = 2
				)
			);
			if ( is( res, "try-error" ) )
			{
				stop( "[MG::PLOT::MULTIPLOT] Error while drawing plot '", plot.type, "': ", geterrmessage() );
			}

			# Updates legend infos
			ltypes <- append( ltypes, ltype );
			plot.distrs <- append( plot.distrs, distr );
		}
		leg.place <- "bottomright";
	} else if ( plot.type == MG_CONSTS_PDF_PLOT ) {
		x <- sort( x );

		# sanity check
		if ( length(x) == 0 )
		{
			warning( "[MG::PLOT::MULTIPLOT] Not enough data to plot '", plot.type, "'" );
			return(invisible(NULL));
		}

		# Plot density histogram
		mg_plot_densHist(
			x,
			title = ifelse( is.null(title), "PDF", title ),
			xlabel = ifelse( is.null(xlabel), "x", xlabel ),
			ylabel = ifelse( is.null(ylabel), "PDF(x)", ylabel ),
			density.col = colors[ col.idx ]
		);

		ltype <- 1; # line type progess
		for ( distr in names( distrs ) )
		{
			# sanity check
			if ( is.null( distrs[[ distr ]] ) || length( distrs[[ distr ]] ) == 0 )
			{
				next;
			}

			# Updates line type progress
			ltype <- ltype + 1;
			if ( ( ltype %% 6 ) == 1 )
			{
				# avoid solid line (used for data)
				ltype <- ltype + 1;
			}

			# Updates color
			col.idx <- col.idx + 1; 

			# plot CDF
			res <- try(
				lines(
					x,
					do.call( "mg_dists_pdf",  c( list(distr), list(x), distrs[[ distr ]] ) ),
					col = colors[ col.idx ],
					lty = ltype,
					lwd = 2
				)
			);
			if ( is( res, "try-error" ) )
			{
				stop( "[MG::PLOT::MULTIPLOT] Error while drawing plot '", plot.type, "': ", geterrmessage() );
			}

			# Updates legend infos
			ltypes <- append( ltypes, ltype );
			plot.distrs <- append( plot.distrs, distr );
		}
		leg.place <- "topright";
	} else {
		# Fall-back case
		warning( "[MG::PLOT::MULTIPLOT] Unknown plot type '", plot.type, "'" );
		return(invisible(NULL));
	}

	# Draws a grid
	grid();

	# Draws a legend
	legend(
		leg.place,
		c( "data", plot.distrs ),
		lty = ltypes,
		col = colors[ 1:col.idx ]
	);
}
