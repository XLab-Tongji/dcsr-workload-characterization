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

## MG_LRD.R
##
## SUMMARY
##  Long Range Dependence functions.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

## MG_LRD_AGGVAR(x)
##
## SUMMARY
##  Aggregate of Variance method.
##
## REFERENCES
##  * M. S. Taqqu and V. Teverovsky. "On Estimating the Intensity of Long-Range
##    Dependence in Finite and Infinite Variance Time Series" (1998)
##
mg_lrd_aggvar <- function( x, n.levels = floor(log10(n)), min.points = 3, cutoffs = 10^c(0.7, 2.5), ret.m = FALSE, ret.vars = FALSE, ret.log.m = FALSE, ret.log.vars = FALSE )
{
	n <- length(x);

	# The increment to be used in plotting and estimating the slope.  Depends on 
	# the number of levels, and the minimum number of points needed at each level.
	increment <- log10( n / min.points ) / n.levels;

	vars.hat <- c();
	#mm <- 1:floor( log10(n) );
	mm <- floor(10^( 0:n.levels * increment ));
	#mm = mm[mm > 0];
	for ( m in mm ) # we use log10 because plots are done in log10 scale
	{
		k.max <- n %/% m;
		xm <- matrix( x[1:(m*k.max)], nrow = k.max, ncol = m, byrow = TRUE );
		vars.hat <- append( vars.hat, var( rowMeans( xm ) ) );
	}

	lsq.w = trunc( (sign( (mm - cutoffs[1]) * (cutoffs[2] - mm) ) + 1) / 2 ); # wt[i] == 1, if mm[i] >= lsq.cut.off[1] && mm[i] <= lsq.cut.off[2]; 0, otherwise
	lmm <- log10(mm);
	lvv <- log10(vars.hat);
	lsq <- lm( lvv ~ lmm, weights = lsq.w );
	H <- lsq$coefficients[2]/2 + 1; # the Hurst parameter
	names(H) <- NULL;

	res <- list();
	res$H <- H;
	res$log.lm <- lsq;
	res$cutoffs <- cutoffs;
	if ( ret.m )
	{
		res$m <- m;
	}
	if ( ret.vars )
	{
		res$vars <- vars.hat;
	}
	if ( ret.log.m )
	{
		res$log.m <- lmm;
	}
	if ( ret.log.vars )
	{
		res$log.vars <- lvv;
	}
	class(res) <- "mg_lrd_aggvar";

	return( res );
}

## MG_LRD_RS
##
## SUMMARY
##  The R/S method.
##
mg_lrd_rs <- function( x, n.levels = floor(log10(length(x))), min.points = 3, cutoffs = 10^c(0.7, 3.5), ret.m = FALSE, ret.rs = FALSE, ret.log.m = FALSE, ret.log.rs = FALSE )
{
	n <- length(x);
	increment <- log10( n/min.points ) / n.levels;
	mm <- floor( 10^( (0:n.levels)*increment ) ); # lags
	#mm = mm[mm > 1]

	yt <- cumsum(as.numeric(x));
	yt2 <- cumsum(as.numeric(x*x));
	rs <- c();
	for ( m in mm )
	{
		timewin <- 1:m;

		s <- sqrt( yt2[m]/m - (yt[m]/m)^2 );

		z <- yt[timewin] - timewin/m * yt[m];

		rs <- append( rs, (max(z) - min(z))/s );
	}

	lsq.w <- trunc((sign((mm-cutoffs[1])*(cutoffs[2]-mm))+1)/2)
	lmm <- log10(mm);
	lrs <- log10(rs);
	lsq <- lm( lrs ~ lmm, weights = lsq.w );
	H <- lsq$coefficients[2];
	names(H) <- NULL;

	res <- list();
	res$H <- H;
	res$log.lm <- lsq;
	res$cutoffs <- cutoffs;
	if ( ret.m )
	{
		res$m <- mm;
	}
	if ( ret.rs )
	{
		res$rs <- rs;
	}
	if ( ret.log.m )
	{
		res$log.m <- lmm;
	}
	if ( ret.log.rs )
	{
		res$log.rs <- lrs;
	}
	class(res) <- "mg_lrd_rs";

	return( res );
}

## MG_LRD_PERIODOGRAM
##
## SUMMARY
##  Periodograms.
##
mg_lrd_periodogram <- function( x, cutoff = 0.10, ret.freqs = FALSE, ret.pgrams = FALSE, ret.log.freqs = FALSE, ret.log.pgrams = FALSE )
{
	n <- length(x);

	FFT <- Mod( fft(x) )^2 / (2 * pi * n);
	pgrams <- FFT[ 1:(n %/% 2 + 1) ];

	# BEGIN NOTE: use this for data coming from stable processes
        #norm <- sum( (x^2)[1:n] );
	#pgrams <- ( FFT/norm )[ 1:(n %/% 2 + 1)];
	# END NOTE: use this for data coming from stable processes

	ff <- ( pi/n ) * c( 2:( n*cutoff ) ); # frequencies
	pgrams <- pgrams[ 2:( n*cutoff ) ]; # periodograms
	#ff <- pi/n; # frequencies
	#pp <- pgrams; # periodograms
	#lsq.w = c( rep.int( 1, n*cutoff-1 ), rep.int( 0, n-n*cutoff ) );
	lsq.w = NULL;
	lff <- log10(ff);
	lpp <- log10(pgrams);
	lsq <- lm( lpp ~ lff, weights = lsq.w );
	H <- (1 - lsq$coefficients[2])/2;
	names(H) <- NULL;
	ff <- (pi/n)*(1:n);
	pgrams <- FFT;

	res <- list();
	res$H <- H;
	res$cutoff <- cutoff;
	res$log.lm <- lsq;
	if ( ret.freqs )
	{
		res$freqs <- ff;
	}
	if ( ret.pgrams )
	{
		res$pgrams <- pgrams;
	}
	if ( ret.log.freqs )
	{
		res$log.freqs <- log10(ff);
	}
	if ( ret.log.pgrams )
	{
		res$log.pgrams <- log10(pgrams);
	}
	class(res) <- "mg_lrd_periodogram";

	return( res );
}

## MG_LRD_PERIODOGRAM.CUM
##
## SUMMARY
##  Cumulative Periodograms.
##
mg_lrd_periodogram.cum <- function( x, cutoff = 0.10, ret.freqs = FALSE, ret.pgrams = FALSE, ret.log.freqs = FALSE, ret.log.pgrams = FALSE )
{
	n <- length(x);

	FFT <- Mod( fft(x) )^2 / (2 * pi * n);
	pgrams <- FFT[ 1:(n %/% 2 + 1) ];

	# BEGIN NOTE: use this for data coming from stable processes
        #norm <- sum( (x^2)[1:n] );
	#pgrams <- ( FFT/norm )[ 1:(n %/% 2 + 1)];
	# END NOTE: use this for data coming from stable processes

	ff <- ( pi/n ) * c( 1:( (n-1)*cutoff ) ); # frequencies
	pgrams <- cumsum( pgrams[ 2:n ] )[ 1:( (n-1)*cutoff ) ]; # cumulative periodograms
	#lsq.w = c( rep.int( 1, n*cutoff-1 ), rep.int( 0, n-n*cutoff ) );
	lsq.w = NULL;
	lff <- log10(ff);
	lpp <- log10(pgrams);
	lsq <- lm( lpp ~ lff, weights = lsq.w );
	H <- (1 - lsq$coefficients[2])/2;
	names(H) <- NULL;
	ff <- (pi/n)*(1:n);
	pgrams <- cumsum(FFT);

	res <- list();
	res$H <- H;
	res$cutoff <- cutoff;
	res$log.lm <- lsq;
	if ( ret.freqs )
	{
		res$freqs <- ff;
	}
	if ( ret.pgrams )
	{
		res$pgrams <- pgrams;
	}
	if ( ret.log.freqs )
	{
		res$log.freqs <- log10(ff);
	}
	if ( ret.log.pgrams )
	{
		res$log.pgrams <- log10(pgrams);
	}
	#class(res) <- c( "mg_lrd_periodogram.cum", "mg_lrd_periodogram" );
	class(res) <- c( "mg_lrd_periodogram" );

	return( res );
}

