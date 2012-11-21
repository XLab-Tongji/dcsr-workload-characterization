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

## MG_FIT
##
## SUMMARY
##  A collection of fitting methods for several distribution.
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

source( "lib/mg/mg_consts.R" );

mg_fit_methodsNameForDistrs <- function( distrs )
{
	methods <- list();

	for ( d in distrs )
	{
		if ( d == MG_CONSTS_BETA_DIST )
		{
			methods[[ d ]] <- c( MG_CONSTS_MLE_FIT, MG_CONSTS_METHODOFMOMENTS_FIT );
		} else if ( d == MG_CONSTS_BINOMIAL_DIST ) {
			methods[[ d ]] <- MG_CONSTS_MLE_FIT;
		} else if ( d == MG_CONSTS_CAUCHY_DIST ) {
			methods[[ d ]] <- MG_CONSTS_MLE_FIT;
		} else if ( d == MG_CONSTS_CHISQUARED_DIST ) {
			methods[[ d ]] <- MG_CONSTS_MLE_FIT;
		} else if ( d == MG_CONSTS_EXPONENTIAL_DIST ) {
			methods[[ d ]] <- MG_CONSTS_MLE_FIT;
		} else if ( d == MG_CONSTS_FISHERF_DIST ) {
			methods[[ d ]] <- MG_CONSTS_MLE_FIT;
		} else if ( d == MG_CONSTS_FRECHET_DIST ) {
			methods[[ d ]] <- MG_CONSTS_MLE_FIT;
		} else if ( d == MG_CONSTS_GAMMA_DIST ) {
			methods[[ d ]] <- c( MG_CONSTS_MLE_FIT, MG_CONSTS_METHODOFMOMENTS_FIT );
		} else if ( d == MG_CONSTS_GENERALIZEDEXTREMEVALUE_DIST ) {
			methods[[ d ]] <- MG_CONSTS_MLE_FIT;
		} else if ( d == MG_CONSTS_GENERALIZEDPARETO_DIST ) {
			methods[[ d ]] <- MG_CONSTS_MLE_FIT;
		} else if ( d == MG_CONSTS_GUMBEL_DIST ) {
			methods[[ d ]] <- c( MG_CONSTS_MLE_FIT, MG_CONSTS_METHODOFMOMENTS_FIT );
		} else if ( d == MG_CONSTS_HYPEREXPONENTIAL_DIST ) {
			methods[[ d ]] <- MG_CONSTS_EM_FIT;
		} else if ( d == MG_CONSTS_LOGISTIC_DIST ) {
			methods[[ d ]] <- MG_CONSTS_MLE_FIT;
		} else if ( d == MG_CONSTS_LOGNORMAL_DIST ) {
			methods[[ d ]] <- MG_CONSTS_MLE_FIT;
		} else if ( d == MG_CONSTS_MMPP_DIST ) {
			# n/a
			#methods[[ d ]] <- MG_CONSTS_EM_FIT;
		} else if ( d == MG_CONSTS_NORMAL_DIST ) {
			methods[[ d ]] <- MG_CONSTS_MLE_FIT;
		} else if ( d == MG_CONSTS_PARETO_DIST ) {
			methods[[ d ]] <- c( MG_CONSTS_MLE_FIT, MG_CONSTS_HILL_FIT, MG_CONSTS_YULE_FIT );
		} else if ( d == MG_CONSTS_CONTINUOUSPHASETYPE_DIST ) {
			methods[[ d ]] <- MG_CONSTS_METHODOFMOMENTS_FIT;
		} else if ( d == MG_CONSTS_POISSON_DIST ) {
			methods[[ d ]] <- MG_CONSTS_MLE_FIT;
		} else if ( d == MG_CONSTS_REVERSEWEIBULL_DIST ) {
			# n/a
		} else if ( d == MG_CONSTS_SKEWEDSTABLE_DIST
			|| d == MG_CONSTS_STABLE_DIST
		) {
			methods[[ d ]] <- c( MG_CONSTS_MLE_FIT, MG_CONSTS_QUANTILEMETHOD_FIT );
		} else if ( d == MG_CONSTS_STUDENTT_DIST ) {
			methods[[ d ]] <- MG_CONSTS_MLE_FIT;
		} else if ( d == MG_CONSTS_SYMSTABLE_DIST ) {
			methods[[ d ]] <- c( MG_CONSTS_MLE_FIT, MG_CONSTS_QUANTILEMETHOD_FIT );
		} else if ( d == MG_CONSTS_UNIFORM_DIST ) {
			methods[[ d ]] <- MG_CONSTS_MLE_FIT;
		} else if ( d == MG_CONSTS_WEIBULL_DIST ) {
			methods[[ d ]] <- c( MG_CONSTS_MLE_FIT, MG_CONSTS_LEASTSQUARED_FIT );
		}
	}

	return( methods );
}

mg_fit_method <- function( distr, methodName, x, ... )
{
	f <- NULL;
	if ( distr == MG_CONSTS_BETA_DIST )
	{
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			f <- "mg_fit_mleBeta";
		} else if ( methodName == MG_CONSTS_METHODOFMOMENTS_FIT ) {
			f <- "mg_fit_momBeta";
		}
	} else if ( distr == MG_CONSTS_BINOMIAL_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			f <- "mg_fit_mleBinomial";
		}
	} else if ( distr == MG_CONSTS_CAUCHY_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			f <- "mg_fit_mleCauchy";
		}
	} else if ( distr == MG_CONSTS_CHISQUARED_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			f <- "mg_fit_mleChiSquared";
		}
	} else if ( distr == MG_CONSTS_CONTINUOUSPHASETYPE_DIST ) {
		if ( methodName == MG_CONSTS_METHODOFMOMENTS_FIT )
		{
			f <- "mg_fit_momCPH";
		}
	} else if ( distr == MG_CONSTS_EXPONENTIAL_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			f <- "mg_fit_mleExponential";
		}
	} else if ( distr == MG_CONSTS_FISHERF_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			f <- "mg_fit_mleFisherF";
		}
	} else if ( distr == MG_CONSTS_FRECHET_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			#f <- "mg_fit_mleFrechet";
			f <- "mg_dists_frechet.fit.mle";
		}
	} else if ( distr == MG_CONSTS_GAMMA_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			f <- "mg_fit_mleGamma";
		} else if ( methodName == MG_CONSTS_METHODOFMOMENTS_FIT ) {
			f <- "mg_fit_momGamma";
		}
	} else if ( distr == MG_CONSTS_GENERALIZEDEXTREMEVALUE_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			f <- "mg_fit_mleGEV";
		}
	} else if ( distr == MG_CONSTS_GENERALIZEDPARETO_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			#f <- "mg_fit_mleGPD";
			f <- "mg_dists_gpd.fit.mle";
		}
	} else if ( distr == MG_CONSTS_GUMBEL_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			#f <- "mg_fit_mleGumbel";
			f <- "mg_dists_gumbel.fit.mle";
		}
	} else if ( distr == MG_CONSTS_HYPEREXPONENTIAL_DIST ) {
		if ( methodName == MG_CONSTS_EM_FIT )
		{
			f <- "mg_fit_emHyperExp";
		}
	} else if ( distr == MG_CONSTS_LOGISTIC_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			f <- "mg_fit_mleLogistic";
		}
	} else if ( distr == MG_CONSTS_LOGNORMAL_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			f <- "mg_fit_mleLogNormal";
		}
	} else if ( distr == MG_CONSTS_MMPP_DIST ) {
		# n/a
		#methods[[ d ]] <- MG_CONSTS_EM_FIT;
	} else if ( distr == MG_CONSTS_NORMAL_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			f <- "mg_fit_mleNormal";
		}
	} else if ( distr == MG_CONSTS_PARETO_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			f <- "mg_fit_mlePareto";
		} else if ( methodName == MG_CONSTS_HILL_FIT ) {
			f <- "mg_fit_hillPareto";
		} else if ( methodName == MG_CONSTS_YULE_FIT ) {
			f <- "mg_fit_yulePareto";
		}
	} else if ( distr == MG_CONSTS_POISSON_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			f <- "mg_fit_mlePoisson";
		}
	} else if ( distr == MG_CONSTS_REVERSEWEIBULL_DIST ) {
		# n/a
	} else if ( distr == MG_CONSTS_SKEWEDSTABLE_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			f <- "mg_fit_mleSkewedStable";
		} else if ( methodName == MG_CONSTS_QUANTILEMETHOD_FIT ) {
			f <- "mg_fit_quantileSkewedStable";
		}
	} else if ( distr == MG_CONSTS_STABLE_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			#f <- "mg_fit_mleStable";
			f <- "mg_dists_stable.fit.mle";
		} else if ( methodName == MG_CONSTS_QUANTILEMETHOD_FIT ) {
			#f <- "mg_fit_quantileStable";
			f <- "mg_dists_stable.fit.quantile";
		}
	} else if ( distr == MG_CONSTS_STUDENTT_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			f <- "mg_fit_mleStudentT";
		}
	} else if ( distr == MG_CONSTS_SYMSTABLE_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			f <- "mg_fit_mleSymStable";
		} else if ( methodName == MG_CONSTS_QUANTILEMETHOD_FIT ) {
			f <- "mg_fit_quantileSymStable";
		}
	} else if ( distr == MG_CONSTS_UNIFORM_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			f <- "mg_fit_mleUniform";
		}
	} else if ( distr == MG_CONSTS_WEIBULL_DIST ) {
		if ( methodName == MG_CONSTS_MLE_FIT )
		{
			f <- "mg_fit_mleWeibull";
		} else if ( methodName == MG_CONSTS_LEASTSQUARED_FIT ) {
			f <- "mg_fit_lsqWeibull";
		}
	}

	if ( is.null( f ) )
	{
		stop( c( "[MG::FIT::UTILS::METHOD] Fitting method '", methodName, "' not available for distribution '", distr, "'." ) );
	}

	return( do.call( f, c( list( x ), list( ... ) ) ) );
}
