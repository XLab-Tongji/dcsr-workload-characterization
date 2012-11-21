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

## MG_GOF_QQ
##
## SUMMARY
##  QQ-plot based GoF tests.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

source( "lib/mg/mg_dists.R" );

## MG_GOF_QQ.TEST(x,dist,...)
##
## SUMMARY
##  One-sample GoF test based on QQ-plot.
##
## PARAMS
##  * x: data values.
##  * dist: probability distribution name.
##  * ...: parameters to probabiity distribution function.
##
## RETURN
##  A list with the following fields:
##  * slope: the slope of reference line in the QQ-plot.
##  * int: the intercept of reference line in the QQ-plot.
##  * qqx: x values of QQ-plot. 
##  * qqy: y values of QQ-plot. 
##
mg_gof_qq.test <- function( x, dist, ... )
{
	return( mg_gof_qq.stat( x, dist, ... ) );

}

## MG_GOF_QQ.STAT(x,dist,...)
##
## SUMMARY
##  Returns One-sample Quantile-Quantile statistics.
##
mg_gof_qq.stat <- function( x, dist, ... )
{
        return(
		mg_gof_qq.stat2( 
                        mg_dists_invcdf( dist, ppoints(x), ... )[order(order(x))], # order(oder(x)) == rank(x)
                        x
                )
        );
}

## MG_GOF_QQ.STAT2(x,y)
##
## SUMMARY
##  Returns Two-sample Quantile-Quantile statistics.
##
mg_gof_qq.stat2 <- function( x, y )
{
	x <- sort( x[!is.na(x)] );
	y <- sort( y[!is.na(y)] );
        qx <- quantile(x, c(0.25, 0.75)); # Q1,Q3 for x
        qy <- quantile(y, c(0.25, 0.75)); # Q1,Q3 for y
        slope <- diff(qy)/diff(qx); # dy/dx
        int <- qy[1] - slope * qx[1]; # intercept
	##@{ Experimental
	# Don't use ... if normality assumption don't hold this is useless
	#l.m <- lm( x ~ y ); # linear regression
	##@} Experimental
	##@{ Experimental
	cor.pearson <- cor.test(x,y,method="pearson");
	#cor.kendall <- cor.test(x,y,method="kendall"); #No: also detect non-linear relations
	#cor.spearman <- cor.test(x,y,method="spearman"); #No: also detect non-linear relations
	##@} Experimental
	##@{ Experimental
	#y.ref <- slope*x + int;
	#abs.diff <- sum( abs( y - y.ref ) ); # Absolute area difference
	#rel.diff <- abs.diff / sum(y.ref); # Relative area difference
	y.ref <- slope*x + int;
	n.x <- length(x);
	n.y <- length(y);
	n.y.ref <- n.x;
	area.qq <- abs( x*y - x[1:(n.x-1)]*y[1:(n.y-1)] );
	area.ref <- abs( x*y.ref - x[1:(n.x-1)]*y.ref[1:(n.y.ref-1)] );
	abs.diff <- sum( abs( area.ref - area.qq ) ); # Absolute area difference
	rel.diff <- abs.diff / sum( area.ref ); # Relative area difference
	##@} Experimental

	return(
		list(
			slope = slope, #slope
			int = int, # intercept
			qqx = x,
			qqy = y,
			q1 = c( x = qx[1], y = qy[1] ), # first quartile
			q3 = c( x = qx[2], y = qy[2] ), # third quartile
			#lm = l.m, # linear regression
			cor.pearson = cor.pearson, # Pearson r correlation coefficient
			#cor.kendall = cor.kendall, # Kendall's tau correlation coefficient
			#cor.spearman = cor.spearman # Spearman correlation coefficient
			abs.area.diff = abs.diff,
			rel.area.diff = rel.diff
		)
	);
}
