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

## MG_GOF_PP
##
## SUMMARY
##  PP-plot based GoF tests.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

source( "lib/mg/mg_dists.R" );

## MG_GOF_PP.TEST1(x,dist,...)
##
## SUMMARY
##  One-sample GoF test based on PP-plot.
##
## PARAMS
##  * x: data values.
##  * dist: probability distribution name.
##  * ...: parameters to probabiity distribution function.
##
## RETURN
##  A list with the following fields:
##  * slope: the slope of reference line in the PP-plot.
##  * int: the intercept of reference line in the PP-plot.
##  * ppx: x values of PP-plot. 
##  * ppy: y values of PP-plot. 
##
mg_gof_pp.test <- function( x, dist, ... )
{
        return( mg_gof_pp.stat( x, dist, ... ) );
}

## MG_GOF_QQ.STAT2(x,y)
##
## SUMMARY
##  Returns Two-samples Quantile-Quantile statistics.
##
mg_gof_pp.stat <- function( x, dist, ... )
{
	return(
		mg_gof_pp.stat2(
			mg_dists_cdf( dist, x, ... ),
			ppoints(x)[order(order(x))]
		)
	);
}

## MG_GOF_QQ.STAT2(x,y)
##
## SUMMARY
##  Returns Two-samples Quantile-Quantile statistics.
##
mg_gof_pp.stat2 <- function( x, y )
{
        x <- sort( x[!is.na(x)] );
        y <- sort( y[!is.na(y)] );

        stats <- list();
        stats$p1 <- c(.25, .25);
        stats$p3 <- c(.75, .75);
        #data.qq <- quantile( xx, probs = c( p1, p3 ) );
        #distr.qq <- mg_dists_invcdf( dist, xx, ... );
        #stats$q1 = c( distr.qq[1], data.qq[1] );
        #stats$q3 = c( distr.qq[2], data.qq[2] );
        ## Draw a line passing through 25% and 75% percentile.
        ##   slope = (y2-y1)/(x2-x1)
        ##   int   = y1 - x1*(y2-y1)/(x2-x1)
        #stats$slope <- (data.qq[2] - data.qq[1]) / (distr.qq[2]-distr.qq[1]);
        #stats$int <- data.qq[1] - distr.qq[1]*stats$slope;
	stats$slope = 1;
	stats$int = 0;

        ##@{ Experimental
        stats$cor.pearson <- cor.test(x,y,method="pearson");
        #stats$cor.kendall <- cor.test(x,y,method="kendall"); #No: also detect non-linear relations
        #stats$cor.spearman <- cor.test(x,y,method="spearman"); #No: also detect non-linear relations
        ##@} Experimental
        ##@{ Experimental
	y.ref <- stats$slope*x + stats$int;
	n.x <- length(x);
	n.y <- length(y);
	n.y.ref <- n.x;
	area.pp <- abs( x*y - x[1:(n.x-1)]*y[1:(n.y-1)] );
	area.ref <- abs( x*y.ref - x[1:(n.x-1)]*y.ref[1:(n.y.ref-1)] );
	stat$abs.area.diff <- sum( abs( area.ref - area.pp ) ); # Absolute area difference
	stats$rel.area.diff <- stats$abs.area.diff / sum( area.ref ); # Relative area difference
        ##@} Experimental

        return( stats );
}
