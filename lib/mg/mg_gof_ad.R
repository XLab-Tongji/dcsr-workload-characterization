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

## MG_GOF_AD
##
## SUMMARY
##  The Anderson-Darling Goodness-of-Fit test.
##
## DESCRIPTION
##  The Anderson-Darling Goodness-of-Fit test.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##
## REFERENCES
##  [1] Anderson T. W. and Darling D. A., "Asymptotic Theory of Certain
##      'Goodness of Fit' Criteria Based on Stochastic Processes", 1952
##  [2] Marsaglia G., "Evaluating the Anderson-Darling Distribution", Journal
##      of Statistical Software, Vol. 9, Issue 2, 2004.
##  [3] Stephens M. A., "EDF Statistics for Goodness of Fit and Some
##      Comparisons", Journal of the American Statistical Association,
##      Vol. 69, Issue 347, 1974.
##

source( "lib/mg/mg_consts.R" );
source( "lib/mg/mg_dists.R" );
source( "lib/mg/mg_eda.R" );
source( "lib/mg/mg_math_poly.R" );
source( "lib/mg/mg_matrix.R" );
source( "lib/mg/mg_utils.R" );

## MG_GOF_AD.TEST(x,dist,...,alternative=c("two.sided", "less", "greater"))
##
## SUMMARY
##  The One-Sample Anderson-Darling test.
##
## PARAMETERS
##  * x: a numeric vector of data values.
##  * dist: the name of theoretical distribution (see 'mg_consts.R').
##  * ...: parameters to pass to CDF function related to the theoretical
##         distribution (e.g. the mean and the standard error in case of Normal
##         distribution). If not specified, the parameters will be estimated
##         from the given data values 'x' (e.g. by means MLE method).
##  * alternative: type of the alternative hyphotesis.
##
## RETURN
##  A list with folling field:
##  * statistc: the test statistic A^2
##  * statistic: the test statistc A^2,
##  * statistic.adj: the adjustement done to the test statistic A^2,
##  * p.value: the p-value (i.e. the minimum significance level above that the
##    null hypothesis has to be refused. For significance levels less than the
##    p-value, the null hypothesis cannot be refused).
##  * p.value.isexact: TRUE if p-value is an exact value.
## 
## DESCRIPTION
##  Performs the Anderson-Darling goodness-of-fit test against the given
##  samples 'x' and distribution 'dist'.
##  The Anderson-Darling test is defined as:
##   A^2 = -N - S
##  where:
##   * N: the sample size
##   * S: \sum_{k=1}^N \frac{2*k-1}{N} \left[\ln F(x_k) + \ln(1 - F(x_{N+1-k})) \right]
##  See [2,3].
##
mg_gof_ad.test <- function( x, dist, ..., alternative = c("two.sided", "less", "greater") )
{
	# >>> TODO: at the moment only the two.sided case is handled <<<

print( "BEFORE [AD]" );
print( list(...) );
print( "AFTER [AD]" );
	alternative <- match.arg( alternative );

	## Two-Sided Anderson-Darling

	x <- sort( na.omit(x) );
	n <- length( x );
	useADCdf <- FALSE;
	adj <- 1; 

	# Compute adjustment and critical values to use for stats.
	if ( dist == MG_CONSTS_NORMAL_DIST )
	{
		# taken from [Stephen1974EDF]

		# Note that the values from NIST dataplot don't work nearly as well.
		adj <- 1 + ( 0.75 + 2.25 / n ) / n;
		#pvals <- c( 0.1, 0.05, 0.025, 0.01 );
		pvals <- c( 0.15, 0.1, 0.05, 0.025, 0.01 );
		Acritn <- matrix(
			c(
				0, .514, .578, .683, .779, .926, # n <= 10
				11, .528, .591, .704, .815, .969, # 10 < n <= 20
				21, .546, .616, .735, .861, 1.021, # 20 < n <= 50
				51, .559, .631, .754, .884, 1.047, # 50 < n <= 100
				101, .576, .656, .787, .918, 1.092 # n > 100
			),
			nrow = 5,
			ncol = 6,
			byrow = TRUE
		);
		#Acrit <- c( 0.631, 0.752, 0.873, 1.035 );
		Acrit <- Acritn[ mg_utils_lookup( Acritn[,1], n, leftmost.closed = TRUE, rightmost.closed = TRUE ) , 2:6 ];

		x <- pnorm( x, ... );
	} else if ( dist == MG_CONSTS_UNIFORM_DIST ) {
		# Put invalid data at the limits of the distribution
		# This will drive the statistic to infinity.
		params <- list(...);
		if ( length(params) )
		{
			# normalize between [0,1]

			min <- params$min;
			max <- params$max;
			x <- (x - min)/(max - min);
		}
		x[ x < 0 ] <- 0; # sanity check
		x[ x > 1 ] <- 1; # sanity check
		pvals <- c( 0.1, 0.05, 0.025, 0.01 );
		Acrit <- c( 1.933, 2.492, 3.070, 3.857 );
		useADCdf <- TRUE;
	} else if ( dist == MG_CONSTS_WEIBULL_DIST && !missing( ... ) ) {
		# taken from http://www.statisticalengineering.com/goodness.htm

		adj <- 1 + 0.2/sqrt(n);
		pvals <- c( 0.1, 0.05, 0.025, 0.01 );
		Acrit <- c( 0.637, 0.757, 0.877, 1.038 );

		useADCdf = TRUE;

		x <- pweibull( x, ... );
	} else if ( dist == MG_CONSTS_EXPONENTIAL_DIST ) {
		# taken from [Stephens1974EDF]

		adj <- 1 + 0.6 / n;
		pvals <- c( 0.15, 0.1, 0.05, 0.025, 0.01 );
		# Critical values depend on n.  Choose the appropriate critical set.
		# These values come from NIST dataplot/src/dp8.f.
		Acritn <- matrix(
			c(
				0, 0.887, 1.022, 1.265, 1.515, 1.888, # n <= 10
				11, 0.898, 1.045, 1.300, 1.556, 1.927, # 10 < n <= 20
				21, 0.911, 1.062, 1.323, 1.582, 1.945, # 20 < n <= 50
				51, 0.916, 1.070, 1.330, 1.595, 1.951, # 50 < n <= 100
				101, 0.922, 1.078, 1.341, 1.606, 1.957 # n > 100
			),
			nrow = 5,
			ncol = 6,
			byrow = TRUE
		);
		Acrit <- Acritn[ mg_utils_lookup( Acritn[,1], n, leftmost.closed = TRUE, rightmost.closed = TRUE ) , 2:6 ];

		x <- pexp( x, ... );
	} else {
		# Fallback to Uniform(0,1) case, i.e. calculate to CDF of
		# the theoretical distribution and consider it as U(0,1)
		# variate.

		# Put invalid data at the limits of the distribution
		# This will drive the statistic to infinity.
		x <- mg_dists_cdf( dist, x, ... );
		x[ x < 0 ] <- 0; # sanity check
		x[ x > 1 ] <- 1; # sanity check
		pvals <- c( 0.1, 0.05, 0.025, 0.01 );
		Acrit <- c( 1.933, 2.492, 3.070, 3.857 );

		useADCdf <- TRUE;
	}

	A2 <- mg_gof_ad.stat( x );
	AA2 <- A2 * adj;

	pval <- 0;

	if ( useADCdf )
	{
		pval <- 1 - mg_gof_ad.cdf( AA2, n );
#	} else if ( dist == MG_CONSTS_NORMAL_DIST ) {
#		# taken from "nortest" R package.
#		if (AA2 < 0.2)
#		{
#			pval <- 1 - exp(-13.436 + 101.14 * AA2 - 223.73 * AA2^2)
#		} else if (AA2 < 0.34) {
#			pval <- 1 - exp(-8.318 + 42.796 * AA2 - 59.938 * AA2^2)
#		} else if (AA2 < 0.6) {
#			pval <- exp(0.9177 - 4.279 * AA2 - 1.38 * AA2^2)
#		} else {
#			pval <- exp(1.2937 - 5.709 * AA2 + 0.0186 * AA2^2)
#		}
#	} else if ( dist == MG_CONSTS_WEIBULL_DIST ) {
##FIXME: don't work
##		# taken from Nelson L. S., "The Anderson-Darling Test for Normality", Journal of Quality Technology, 1998
##		pval <- 3.6789468*exp(-AA2/0.1749916);
#
##FIXME: don't work
##		# taken from Romen J. L., "Anderson-Darling: A Goodness of Fit Test for Small Samples Assumptions", START, 2003
##		pval <- 1/(1+exp(-0.1+(1.24*log(AA2))+(4.48*AA2)));
#
#		# fallback to default case (use this while one of the above alternative is fixed)
#		#idx <- mg_utils_lookup( c( -Inf, Acrit ), AA2, leftmost.closed = TRUE, rightmost.closed = TRUE );
#		#pval <- c( 1, pvals )[ idx ];
#
#		# interpolates critical values
#		pval <- approx( Acrit, y = pvals, xout = AA2, yleft = 1, yright = 0 );
	} else {
		# Use critical values table "Acrit"

		# OLD METHOD: bin critical values and get the one falling in the right interval
		#idx <- mg_utils_lookup( c( -Inf, Acrit ), AA2, leftmost.closed = TRUE, rightmost.closed = TRUE );
		#pval <- c( 1, pvals )[ idx ];

		# interpolates critical values
		pval <- approx( Acrit, y = pvals, xout = AA2, yleft = 1, yright = 0 )[["y"]];
	}

	return(
		list(
			statistic = A2,
			statistic.adj = adj,
			p.value = pval,
			p.value.isexact = isTRUE( useADCdf ),
			alternative = alternative
			#critical.values = matrix( c( 100 * (1 - pvals), Acrit ), nrow = length(pvals), byrow = FALSE )
		)
	);
}

## MG_GOF_AD.TEST.BOOT.PAR(x,dist,...,alternative=c("two.sided", "less", "greater"))
##
## SUMMARY
##  The One-Sample Anderson-Darling test with Parametric Bootstraping.
##
## PARAMETERS
##  * x: a numeric vector of data values.
##  * dist: the name of theoretical distribution (see 'mg_consts.R').
##  * ...: parameters to pass to CDF function related to the theoretical
##         distribution (e.g. the mean and the standard error in case of Normal
##         distribution). If not specified, the parameters will be estimated
##         from the given data values 'x' (e.g. by means MLE method).
##  * n.boot: number of bootstrap iterations.
##  * alternative: type of the alternative hyphotesis.
##
## RETURN
##  A list with folling field:
##  * statistc: the test statistic A^2
##  * statistic: the test statistc A^2,
##  * statistic.adj: the adjustement done to the test statistic A^2,
##  * p.value: the p-value (i.e. the minimum significance level above that the
##    null hypothesis has to be refused. For significance levels less than the
##    p-value, the null hypothesis cannot be refused).
##  * p.value.isexact: TRUE if p-value is an exact value.
##  * p.value.boot: p-value obtained from bootstrap method.
##  * alternative: type of alternative hypothesis
##  * statistic.boot.mean: mean of A^2 bootstrap statistic.
##  * statistic.boot.sd: standard deviation of A^2 bootstrap statistic.
## 
## DESCRIPTION
##  Performs the Anderson-Darling goodness-of-fit test against the given
##  samples 'x' and distribution 'dist'.
##  The Anderson-Darling test is defined as:
##   A^2 = -N - S
##  where:
##   * N: the sample size
##   * S: \sum_{k=1}^N \frac{2*k-1}{N} \left[\ln F(x_k) + \ln(1 - F(x_{N+1-k})) \right]
##  See [2,3].
##
mg_gof_ad.test.boot.par <- function( x, dist, ..., n.boot = 999, alternative = c("two.sided", "less", "greater"), fit.methods = NULL )
{
        n.x <- length(x);

        # Estimates parameters for dist from X^* with the same method used
        # for parameters passed as argument
        if ( is.null( fit.methods ) || length( fit.methods ) == 0 )
        {
                fit.methods <- mg_fit_methodsNameForDistrs( dist );
        } else {
                tmp <- fit.methods;
                fit.methods <- list();
                fit.methods[[ dist ]] <- tmp;
        }

        # Perform a A-D test from X and distr (whose parameters are estimated from X)
        ad <- try( mg_gof_ad.test( x, dist, ..., alternative = alternative ) );
	if ( is( ad, "try-error" ) )
	{
		stop( "[MG::GOF::AD::TEST::BOOT::PAR] Error while performing A-D test for distribution '", dist, "' on original sample: ", geterrmessage() );
	}

        ad.pval <- 0;
        boot.mean <- 0;
        boot.var <- 0;
        for ( i in 1:n.boot )
        {
                # Samples from the theoretical distribution (=> the bootstrap sample X^*)
                x.star <- try( mg_dists_rvg( dist, n.x, ... ) );
                if ( is( x.star, "try-error" ) )
                {
                        warning( "[MG::GOF::AD.TEST.BOOT.PAR] Error while sampling from distribution '", dist, "' in bootstrap iteration #", i, ": ", geterrmessage(), "." );
                        next;
                }

                theta.star <- try( mg_fit_multifit( x.star, fit.methods ) );
                if ( is( theta.star, "try-error" ) )
                {
                        warning( "[MG::GOF::AD.TEST.BOOT.PAR] Error while fitting distribution '", dist, "' in bootstrap iteration #", i, ": ", geterrmessage(), "." );
                        next;
                } else if ( is.null( theta.star ) || length( theta.star ) == 0 ) {
                        warning( "[MG::GOF::AD.TEST.BOOT.PAR] Error while fitting distribution '", dist, "' in bootstrap iteration #", i, ": returned NULL." );
                        next;
		}

                # Performs a A-D test between the bootstrap sample X^* and the distribution dist with parameters theta^*
                ad.star <- try( do.call( "mg_gof_ad.test", c( list(x = x.star, dist = dist), theta.star[[dist]]$estimate , list(alternative = alternative) ) ) );
		if ( is( ad.star, "try-error" ) )
		{
			stop( "[MG::GOF::AD::TEST::BOOT::PAR] Error while performing A-D test for distribution '", dist, "' on ", i, "-th bootstrap sample: ", geterrmessage() );
		}

                boot.mean <- ad.star$statistic;
                boot.var <- ad.star$statistic^2;

                if ( ad.star$statistic >= ad$statistic )
                {
                        ad.pval <- ad.pval + 1;
                }
        }
        ad.pval <- ad.pval / (n.boot + 1);
        boot.mean <- boot.mean/n.boot;
        boot.var <- 1/(n.boot-1) * (boot.var - n.boot*boot.mean^2);

        return(
                list(
                        statistic = ad$statistic,
                        statistic.adj = ad$statistic.adj,
                        p.value = ad$p.value,
                        p.value.isexact = ad$p.value.isexact,
                        p.value.boot = ad.pval,
                        statistic.boot.mean = boot.mean,
                        statistic.boot.sd = sqrt(boot.var),
                        alternative = ad$alternative
                )
        );
}

## MG_GOF_AD.STAT
##
## SUMMARY
##  Returns the Anderson-Darling statistic A_n^2.
##
## PARAMETERS
##  * u: a set of ordered uniform [0,1] variates.
##
## RETURN
##  The value of the Anderson-Darling statistic A_n^2.
##
## DESCRIPTION
##  Calculates and returnes the Anderson-Darling statistic:
##   A_n^2 = - n - \frac{1}{n} \sum_{j=1}^n\{(2j - 1)[\ln(U_j) + \ln(1 - U_j)]\}
##  where U_1<U_2<...<U_n is an ordered set of purported uniform [0,1) variates
##  (e.g. in A-D test, U_i are CDF values).
##
mg_gof_ad.stat <- function( u )
{
	n <- length( u );
	#n <- mg_matrix_size( u, 1 );

	if ( n == 0 )
	{
		return( 0 );
	}

#	S <- 0;
#	for (k in 1:n)
#	{
#		S <- S + (2*k-1)/n * ( log(x[k]) + log(1-x[n+1-k] ) );
#	}
#	A2 <- -n - S;
	S <- ( 2 * seq( from = 1, to = n ) - 1 ) * ( log( u ) + log( 1 - rev( u ) ) );
	A2 <- -n - mean(S);
	#i <- mg_matrix_vtprod( 1:n, mg_matrix_ones( 1, mg_matrix_size( u, 2 ) ) );
	#A2 <- -n - sum( (2 * i - 1) * ( log(u) + log(1 - u[n:1]) ) ) / n;

	return( A2 );
}

## MG_GOF_AD.cdf( A2, n )
##
## SUMMARY
##  Returns the CDF of the Anderson-Darling distribution.
##
## PARAMS
##  * A2: the value of the Anderson-Darling statistic
##  * n: the number of samples.
##
## RETURN
##  The CDF \Pr{A_n^2 <= A2\}.
##
## DESCRIPTION
##  Calculates and returns the CDF of the Anderson-Darling distribution for the
##  Anderson-Darling statistic value @c A2: \Pr\{A_n^2 <= A2\}
##
mg_gof_ad.cdf <- function( A2, n )
{
	y <- mg_gof_ad._Ainf( A2 );
	y <- y + mg_gof_ad._errfix( y, n );

	return( y );
}

## MG_GOF_AD._AINF
##
## SUMMARY
##  Asymptotic distribution of Anderson-Darling (A_n^2) statistic.
##
## PARAMS
##  * z: an ordered set of purported uniform [0,1] variates.
##
## RETURN
##  The A_{\infty}^2 statistic.
##
## DESCRIPTION
##  Evaluates the limiting distribution Ainf of the Anderson-Darling statistic:
##   A_n^2 = - n - \frac{1}{n} \sum_{j=1}^n\{(2j - 1)[\ln(U_j) + \ln(1 - U_j)]\}
##  where U_1<U_2<...<U_n is an ordered set of purported uniform [0,1) variates.
##  The function is A_{\infty}^2 (z)=lim_{n->\infty} Pr\{A_n^2 <= z\}.
##  This function calculates an approximation of the real Ainf value.
##  See [2].
##
mg_gof_ad._Ainf <- function( z )
{
	y <- length( z );
	#y <- mg_matrix_zeros( mg_matrix_size( z ) );

	idx <- which( z < 2 );
	if ( any(idx) )
	{
		p <- c( .00168691, -.0116720, .0347962, -.0649821, .247105, 2.00012 );
		z1 <- z[idx];
		y[idx] <- exp( -1.2337141 / z1 ) / sqrt(z1) * mg_math_polyval(p,z1);
	}

	idx <-which( z >= 2 );
	if ( any(idx) )
	{
		p <- c( -.0003146, +.008056, -.082433, +.43424, -2.30695, 1.0776 );
		y[idx] <- exp( -exp( mg_math_polyval(p, z[idx] ) ) );
	}

	return( y );
}

## MG_GOF_AD._ERRFIX( x, n )
##
## SUMMARY
##  Correction for the asymptotic Anderson-Darling statistic calculated with
##  the @c mg_gof_ad._Ainf function.
##
## DESCRIPTION
##  Corrects the error caused by using the asymptotic approximation,
##  mg_gof_ad._Ainf(z).
##  Thus x+ad._errfix(x,n) is uniform in [0,1) for practical purposes;
##  accuracy may be off at the 5th, rarely at the 4th, digit.
##  See [2].
##
mg_gof_ad._errfix <- function( x, n )
{
	if ( is.numeric(n) )
	{
		n <- rep( n, length(x) );
		#n <- n * mg_matrix_ones( mg_matrix_size(x) );
	} else if ( is.numeric(x) ) {
		x <- rep( x, length(n) );
		#x <- x * mg_matrix_ones( mg_matrix_size(x) );
	}

	ret <- mg_vector_zeros( length(x) );
	#ret <- mg_matrix_zeros( mg_matrix_size(x) );
	c <- .01265 + .1757 / n; # the c(n) function

	idx = which( x >= 0.8 );
	if ( any(idx) )
	{
		p <- c( 255.7844, -1116.360, 1950.646, -1705.091, 745.2337, -130.2137 );
		g3 <- mg_math_polyval( p, x[idx] ); # the g3(x) function
		ret[idx] <- g3 / n[idx];
	}

	idx = which( x < 0.8 & x > c );
	if ( any(idx) )
	{
		p = c( 1.91864, -8.259, 14.458, -14.6538, 6.54034, -.00022633 );
		n1 <- 1 / n[idx];
		c1 <- c[idx];
		g2 <- mg_math_polyval( p, (x[idx] - c1) / (0.8 - c1) ); # the g2(x) function
		ret[idx] <- (.04213 + .01365*n1)*n1 * g2;
	}

	idx <- which( x <= c );
	if ( any(idx) )
	{
		x1 <- x[idx] / c[idx]
		n1 <- 1 / n[idx];
		g1 <- sqrt(x1) * (1 - x1) * (49*x1 - 102); # the g1(x) function
		ret[idx] <- ( (0.0037*n1 + 0.00078) * n1 + 0.00006 ) * n1 * g1;
	}

	return( ret );
}

#mg_gof_ad._testunit <- function()
#{
##	x <- matrix( c( 10*rexp(12*10000) ), nrow=12,ncol=10000 );
##	ad <- mg_gof_ad.test( x, MG_CONSTS_EXPONENTIAL_DIST, 1/mean(x) );
##	c <- ad$p.value;
##	m <- matrix( c( 100*c, 100*c( unique(c),1 ) ), nrow=length(c), byrow=FALSE );
##	print( m );
#	r <-c(
#		3.6029739, 0.1450064, 2.6972021, 0.9234109, 0.5358326, 1.2282242, 0.1917170, 0.2979628, 0.7604708, 0.7365290,
#		0.1386191, 0.3310959, 1.0294556, 1.6095886, 1.3369319, 0.1597857, 1.4035689, 0.1648368, 1.0559846, 0.1759148,
#		0.6743113, 0.9745661, 0.3863615, 0.0098597, 1.9261793, 0.1768004, 0.3091252, 1.4070909, 2.3438880, 1.9321391,
#		0.4512831, 0.5376761, 0.8744115, 0.0197727, 3.5359983, 0.0197481, 2.4517625, 0.7649851, 0.0496341, 0.2780023,
#		0.7255745, 0.7740526, 0.1136753, 4.3205383, 0.2717797, 0.9489695, 1.1390225, 0.0169961, 0.7420047, 1.3729246,
#		0.3415969, 0.0838968, 0.8700157, 0.6995485, 0.6041986, 0.2191400, 0.5594274, 0.1409002, 0.7781627, 1.0041561,
#		0.2415441, 2.4736436, 0.9496428, 0.4940230, 0.5989277, 2.1133386, 0.4783774, 0.8229692, 0.7280467, 1.7735080,
#		0.4643371, 0.6883693, 0.3805753, 3.3046616, 1.1859609, 4.3278137, 1.3363834, 1.5032329, 0.0622893, 0.2675796,
#		1.1493986, 0.4865569, 1.1035029, 5.2019019, 0.6063860, 0.2928286, 0.1086893, 2.7934619, 0.0151364, 0.4575336,
#		0.7923198, 0.8104809, 0.5169460, 0.9808626, 1.7326085, 2.0026444, 0.2407728, 0.2196102, 0.4462588, 1.1612987,
#		0.0576662, 2.8035915, 1.0835594, 0.0201025, 0.9273553, 0.2917325, 0.5558468, 0.6111510, 0.1863934, 0.0184351,
#		0.2407093, 2.1243305, 0.6580299, 0.8583535, 2.0350314, 0.3420737, 0.1286765, 0.0310931, 1.3868716, 0.8793901
#	);
#	r <- matrix( r, nrow=12, ncol=10, byrow=T );
#
#	set.seed(1);
#	cc <- c();
#	for ( k in c(1:10) )
#	{
#		#x <- 10*rexp(12);
#		x <- 10*r[,k];
#		ad <- mg_gof_ad.test( x, MG_CONSTS_EXPONENTIAL_DIST, 1/mean(x) );
#		print( ad );
#		c <- ad$p.value;
#		cc <- append( cc, 100*c );
#		#m <- matrix( c( 100*c, 100*c( unique(c),1 ) ), nrow=length(c), byrow=FALSE );
#		#print( m );
#	}
#	append( cc, 1 );
#	ccf <- factor( cc );
#
#	print( "[Bin ; AbsFreq] table:" );
#	print( table(ccf) );
#
#	# The output should be:
#	# 2.5  10 100 
#	#   1   1   8 
#
#}
