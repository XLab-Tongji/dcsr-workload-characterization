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

## MG_EDA
##
## SUMMARY
##  Exploratory Data Analysis.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

mg_eda_ecdf <- function( x )
{
	Fn <- ecdf( x );
	class(Fn) <- c("mg_eda_ecdf", class(Fn));

	return( Fn );
}

mg_eda_eccdf <- function ( x )  
{
	#x <- sort(x) 
	n <- length(x);
	if ( n < 1 )
	{
		stop("[MG::EDA::ECCDF] 'x' must have 1 or more non-missing values");
	}
	vals <- sort(unique(x));
	rval <- approxfun(
			vals,
			1-cumsum(tabulate(match(x, vals)))/n,
			method = "constant",
			yleft = 1,
			yright = 0,
			f = 0,
			ties = "ordered"
	);
	class(rval) <- c("mg_eda_eccdf", "eccdf", "stepfun", class(rval));
	attr(rval, "call") <- sys.call();

	return( rval );
}

## MG_EDA_SPLIT(x,n.groups)
##
## SUMMARY
##  Splits the vector @c x into @c n.groups groups
##
mg_eda_split <- function( x, n.groups )
{
	x.len <- length(x);
	n.groups <- n.groups + x.len %% n.groups; # adjusted number of groups

	if ( n.groups == 1 )
	{
		# No split
		return( as.matrix(x) );
	}

	group.size <- ceiling( x.len / n.groups );
	groups <- matrix( NA, nrow = n.groups, ncol = group.size );
	for ( i in 1:n.groups )
	{
		begin <- (i-1) * group.size + 1;
		#end <- min( begin + (group.size - 1), x.len );
		end <- begin + (group.size - 1); # NA automatically added

		groups[i,] <- x[ begin:end ];
	}
	return( groups );
}

## MG_EDA_AGGREGATE(x,n.groups,fit)
##
## SUMMARY
##  Aggregates a data set.
##
mg_eda_aggregate <- function( x, n.groups, fit.method = c("none", "accumulate", "enlarge") )
{
	fit.method <- match.arg(fit.method);

	x.len <- length(x);
	if ( x.len == n.groups )
	{
		# No grouping
		return( x );
	}

	have.residuals <- FALSE;
	#n.groups <- n.groups + x.len %% n.groups; # adjusted number of groups
	#group.size <- ceiling( x.len / n.groups );
	group.size <- 0;
	if ( x.len %% n.groups )
	{
		have.residuals <- TRUE;

		# n.groups is not a divisor of x.len
		if ( fit.method == "none" )
		{
			# No fitting
			# ==> residual elements will be discarded

			group.size <- floor( x.len / n.groups );
		} else if ( fit.method == "accumulate" ) {
			# Accumulate fitting
			# ==> The number of groups must remain fixed.
			#     The last group will be bigger than others.
			#     (handled below)

			group.size <- floor( x.len / n.groups );
		} else {
			# Enlarge fitting
			# ==> The number of groups will be enlarged accordingly

			group.size <- ceiling( x.len / n.groups );
			n.groups <- x.len / group.size;
		}
	} else {
		# n.groups is a divisor of x.len

		group.size <- x.len / n.groups;
	}
	for ( i in 1:n.groups )
	{
		begin <- (i-1) * group.size + 1;
		end <- begin + (group.size - 1); # NA automatically added

		x[i] <- sum( x[ begin:end ] );
	}
	if ( have.residuals && fit.method == "accumulate" )
	{
		# Fixed fitting
		# => Last group get all remaining elements.

		x[n.groups] <- x[n.groups] + sum( x[(n.groups*group.size+1):x.len] );
	}
	x <- x[1:n.groups];
	return( x );
}

mg_eda.outliers.grubbs <- function(x, alternative=c("two.sided","max","min") )
{
	haveLib <- require("outliers");
	if ( !haveLib )
	{
		stop( "[MG::EDA::OUTLIIERS::GRUBBS] Missing dependencies." );
	}

	alternative = match.arg(alternative);

	#if ( alternative == "two.sided" )
	#{
	#	type <- 
	noout <- FALSE; # TRUE if there're no outlier
	while ( length(x) > 0 && !noout )
	{
		g <- grubbs.test(x);
		#TODO
	}
	
}

## MG_EDA.SHAPE.PLT
##
## SUMMARY
##  Plots Shape informations.
##
mg_eda.shape.plot <- function(x)
{
	par(mfrow = c(2, 2));
	mg_plot_freqHist( x );
	grid();
	boxplot(x,col="gray90",border="royalblue",pch="+", main="Box Plot");
	grid();
	#iqd <- summary(x)[5] - summary(x)[2];
	iqd <- IQR(x);
	#plot(density(x, width = 2 * iqd), xlab = "x", ylab = "", type = "l",col="royalblue",lwd=2);
	mg_plot_densHist(x, density.width = 2*iqd, title = "Density Histogram");
	grid();
	qqnorm(x,col="royalblue",pch="+");
	qqline(x,col="black",lty=2);
	grid();

	return( invisible() );
}

## MG_EDA.ACOR.PLT
##
## SUMMARY
##  Plots Autocorrelation informations.
##
mg_eda.acor.plot <- function(x)
{
	par(mfrow = c(2, 1));
	ts.plot( as.ts(x), type = "h", col = "royalblue", main = "Run Sequence Plot", xlab = "index", ylab = "x" );
	grid();
	acf( x, main = "ACF Plot" );
	grid();

	return( invisible() );
}

##    F_i = i/n\,
##    S_i = \Sigma_{j=1}^i \; y_j\,
##    L_i = S_i / S_n \, 
mg_eda.lorenz <- function(x)
{
	x <- sort(x);
	n <- length(x);
	f <- (1:n)/n;
	s <- cumsum(x);
	l <- s/s[length(s)];
	#G <- mg_eda.gini(x);
	#plot( c(0,f), c(0,l), type = "l" );

	#return( invisible( list( x = c(0,x), count = c(0,f), mass = c(0,l), gini = G ) ) );
	return( invisible( list( x = c(0,x), count = c(0,f), mass = c(0,l) ) ) );
}

## MG_EDA.GINI(x)
##
## SUMMARY
##  Returns the Gini coefficient
mg_eda.gini <- function(x)
{
#	n <- length(x);
#	x <- sort(x);
#	G <- sum(x * 1:n);
#	G <- 2*G/(n*sum(x));
#	return( G - 1 - (1/n) );

	n <- length(x);
	G <- 0;
	for ( i in 1:n )
	{
		G <- G + sum( abs( x[i] - x ) );
	}
	G <- G / (2 * n * sum(x));

	return( G );
}

#################################################
# percentile computation functions, adapted from
# SAS Procedures Guide, Version 6, 3rd ed., 625f.
#
## for testing:
#test.x <- c (1,2,3,4,5,6,7,8,9,4,2,3,6,7,8,2,5,2)
#test.p <- c (26, 51)
#for (i in 1:5) print (perc (test.x, test.p, i))
## end.
mg_eda.percentile <- function (x, p = c(5,10,25,50,75,95), pctldef = 5)
{
	#
	# initialization lacking error checking
	#
	x <- na.omit (as.vector (x)) # transform input into non-empty vector
	p <- na.omit (as.vector (p)) # wanted percentiles into non-empty vector
	x <- sort    (x)             # sort it in ascending order
	n <- length  (x)             # length of non-empty vector
	#
	if (pctldef == 4) n <- n + 1 # increase n by one in case pctldef == 4
	#
	j <- trunc   (n * p / 100)   # j is the integer part of the product n * p
	g <- - j   + (n * p / 100)   # g is the fractional part of the product n * p
	#
	# the different computational procedures follow
	#
	if (pctldef == 1)            # weighted average at x\sub{np}
	{
		# x\sub{0} is taken to be x\sub{1}, cf. above
		perc <- (1 - g) * x [j] + g * x [j + 1]
	}
	#
	if (pctldef == 2)  {         # observation number closest to np
		i <- ifelse (g == 0.5, ifelse (trunc (j / 2) == j / 2, j, j + 1), trunc  (n * p / 100 + 0.5))
		perc <- x [i]
	}
	#
	if (pctldef == 3)            # empirical distribution function
	{
		perc <- x [ifelse (g == 0, j, j + 1)]
	}
	#
	if (pctldef == 4)            # weighted average aimed at x\sub{p*(n+1)}
	{ # x\sub{n+1} is taken to be x\sub{n}, cf.  above
		perc <- (1 - g) * x [j] + g * x [j + 1]
	}
	#
	if (pctldef == 5)            # empirical distribution function with averaging
	{
		perc <- ifelse (g == 0, (x [j] + x [j + 1]) / 2, x [j + 1])
	}
	#
	names  (perc) <- p
	return (perc)
}
