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

## MG_APP_TERAGRID-OVERALL-IAT-ALL
##
## SUMMARY
##  TeraGrid - Overall interarrival times.
##
## DESCRIPTION
##  Writes:
##  * a data-file containing the overall interarrival times for TERAGRID log;
##  * a stat (ascii) file containing varoius statistics.
##  Plots:
##  * a frequency histogram;
##  * a density histogram;
##  * a time-series plot;
##  * an LLCD plot with several fitting distributions;
##  * a QQ-plot;
##  * an ACF plot;
##  * a Hill estimator plot.
##  * ...
##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

oldwarn <- getOption( "warn" );
options( warn = 1 );

source( "app/mg_app_conf.R" );
source( "lib/mg/mg_debug.R" );
source( "lib/mg/mg_dists.R" );
source( "lib/mg/mg_eda.R" );
source( "lib/mg/mg_fit.R" );
source( "lib/mg/mg_gof.R" );
source( "lib/mg/mg_teragrid.R" );
source( "lib/mg/mg_lrd.R" );
source( "lib/mg/mg_moments.R" );
source( "lib/mg/mg_plot.R" );
source( "lib/mg/mg_statsutils.R" );

set.seed(123456); # used to eventually repeat the same experiment

dbg <- new( "mg_debug", activate = TRUE );

#teragrid <- new( "mg_teragrid", MG_APP_CONF_TERAGRID_LOGFULLNAME );
teragrid <- new(
	"mg_teragrid",
	MG_APP_CONF_TERAGRID_LOGFULLNAME,
	omitSeqE = TRUE,
	omitSgtE = TRUE,
	omitNle0 = TRUE,
	fixAgtS = TRUE
);

mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Reading IAT file" );

iat <- mg_teragrid_overallInterArrTimes( teragrid );

sink( file = paste( MG_APP_CONF_TERAGRID_STATSPATH, "/teragrid-overall-iat-stats-gof.out", sep = "" ) );

# BEGIN Plot fitted distributions and writes estimated parameters to data file.
mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Fitting distributions" );
fitDistrs <- mg_fit_methodsNameForDistrs(
	c(
#		MG_CONSTS_BETA_DIST,
		MG_CONSTS_CAUCHY_DIST,
#		MG_CONSTS_CHISQUARED_DIST,
		MG_CONSTS_EXPONENTIAL_DIST,
#		MG_CONSTS_FISHERF_DIST,
#		#MG_CONSTS_FRECHET_DIST,
		MG_CONSTS_GAMMA_DIST,
		MG_CONSTS_GENERALIZEDEXTREMEVALUE_DIST,
		MG_CONSTS_GENERALIZEDPARETO_DIST,
#		#MG_CONSTS_GUMBEL_DIST,
#		#MG_CONSTS_HYPEREXPONENTIAL_DIST,
		MG_CONSTS_LOGISTIC_DIST,
		MG_CONSTS_LOGNORMAL_DIST,
#		#MG_CONSTS_MMPP_DIST,
		MG_CONSTS_NORMAL_DIST,
		MG_CONSTS_PARETO_DIST,
		MG_CONSTS_CONTINUOUSPHASETYPE_DIST,
#		#MG_CONSTS_REVERSEWEIBULL_DIST,
#		MG_CONSTS_SKEWEDSTABLE_DIST,
#		MG_CONSTS_STUDENTT_DIST,
#		MG_CONSTS_SYMSTABLE_DIST,
		MG_CONSTS_WEIBULL_DIST
	)
);
if ( length( fitDistrs ) == 0 )
{
	warning( "Cannot perform distributions fitting." );
	cat( "(WARNING) No method available for performing distributions fitting.", "\n" );
	quit();
}
#postscript(
#	file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-fit-llcd.eps", sep = "" ),
#	#title = "TERAGRID - Overall Interarrival Times",
#	horizontal = FALSE,
#	onefile = FALSE
#);
mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Distribution multi fitting" );
fit <- try(
	mg_fit_multifit(
		iat,
		fitDistrs
		#plotType = MG_CONSTS_LLCD_PLOT
	)
);
#dev.off();
if ( is( fit, "try-error" ) )
{
	#unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-fit-llcd.eps", sep = "" ) );
	warning( "Error while performing distributions fitting." );
	warning( "Error message is: ", geterrmessage() );
	quit();
}
cat( "## Distributions Fitting:", "\n" );
print( fit );

##@{ Performs GoF tests

cat( "## Goodness-of-Fit tests:", "\n" );
fitDistrs <- names( fit );
for ( fitDistr in fitDistrs )
{
	cat( paste( "GoF tests for distribution: ", fitDistr, sep = "" ), "\n" );
	if ( length( fit[[ fitDistr ]] ) == 0 )
	{
		cat( "(WARNING) Fitting informations not available.", "\n" );
		next;
	}

	# BEGIN GoF tests
#BEGIN don't work
#	res <- try(
#		mg_gof_multigof.test(
#			iat,
#			dist = fitDistr,
#			tests = c(
#				MG_CONSTS_ANDERSONDARLING_GOF,
#				MG_CONSTS_KOLMOGOROVSMIRNOV_GOF,
#				MG_CONSTS_CHISQUARED_GOF
#			),
#			fit[[ fitDistr ]]
#		)
#	);
#END don't work
	mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Doing Goodness-of-Fit tests for '", fitDistr, "' distribution" );
	res <- try(
		do.call(
			"mg_gof_multigof.test",
			c(
				list(iat),
				list(dist = fitDistr),
				fit[[ fitDistr ]]$estimate,
				list(
					tests = c(
#						MG_CONSTS_ANDERSONDARLING_GOF,
#						MG_CONSTS_ANDERSONDARLINGBOOTPAR_GOF,
#						MG_CONSTS_KOLMOGOROVSMIRNOV_GOF,
#						MG_CONSTS_KOLMOGOROVSMIRNOVBOOTPAR_GOF,
#						#MG_CONSTS_KOLMOGOROVSMIRNOVBOOTNONPAR_GOF,
						#MG_CONSTS_KOLMOGOROVSMIRNOVSIM_GOF,
						#MG_CONSTS_KOLMOGOROVSMIRNOV2_GOF,
						#MG_CONSTS_KOLMOGOROVSMIRNOV2BOOT_GOF,
						MG_CONSTS_CHISQUAREDPAR_GOF
#						MG_CONSTS_QQ_GOF
					),
					dist.nparams = length( fit[[ fitDistr ]]$estimate ),
					boot.n = 999,
					boot.fit.methods = fit[[ fitDistr ]]$methods
				)
			)
		)
	);
	if ( !is( res, "try-error" ) )
	{
		print( res );
	} else {
		warning( "Error while doing GoF tests for distribution '", fitDistr, "'." );
		warning( "Error message is: ", geterrmessage() );
		cat( "(ERROR) GoF tests failed.", "\n" );
	}
	# END GoF tests

	# BEGIN Q-Q Plot
	mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Drawing Q-Q Plot for '", fitDistr, "' distribution" );
	cat( "## Q-Q Plot:", "\n" );
	postscript(
		file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-", fitDistr, "fit-qqplot.eps", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		do.call(
			"mg_plot_qq",
			c(
				list(iat),
				list(dist = fitDistr),
				fit[[ fitDistr ]]$estimate
			)
		)
	);
	dev.off();
	if ( !is( res, "try-error" ) )
	{
		cat( "\t", paste( "Slope: ", res$slope, " - Intercept: ", res$int ), "\n" );
		cat( "\t", "Correlation Coeff.: ", paste( "Pearson: ", res$cor.pearson$estimate, " @ [", res$cor.pearson$conf.int[1], ",", res$cor.pearson$conf.int[2], "] (p-value: ", res$cor.pearson$p.value, ")", sep = "" ), "\n" );
		cat( "\t", "Area Difference: Absolute: ", res$abs.area.diff, ", Relative: ", res$rel.area.diff, "\n" );
	} else {
		unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-", fitDistr, "fit-qqplot.eps", sep = "" ) );
		warning( "Error while drawing Q-Q Plot distribution: '", fitDistr, "'." );
		warning( "Error message is: ", geterrmessage() );
		cat( "(ERROR) Q-Q Plot failed.", "\n" );
	}
	# END Q-Q Plot

#	# BEGIN Log-Log Q-Q Plot
#	mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Drawing Log-Log Q-Q Plot for '", fitDistr, "' distribution" );
#	cat( "## Log-Log Q-Q Plot:", "\n" );
#	postscript(
#		file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-", fitDistr, "fit-llqqplot.eps", sep = "" ),
#		horizontal = FALSE,
#		onefile = FALSE
#	);
#	res <- try(
#		do.call(
#			"mg_plot_qq",
#			c(
#				list(iat),
#				list(dist = fitDistr),
#				fit[[ fitDistr ]]$estimate,
#				log = "xy"
#			)
#		)
#	);
#	dev.off();
#	if ( !is( res, "try-error" ) )
#	{
#		cat( "\t", paste( "Slope: ", res$slope, " - Intercept: ", res$int ), "\n" );
#		cat( "\t", "Correlation Coeff.: ", paste( "Pearson: ", res$cor.pearson$estimate, " @ [", res$cor.pearson$conf.int[1], ",", res$cor.pearson$conf.int[2], "] (p-value: ", res$cor.pearson$p.value, ")", sep = "" ), "\n" );
#	} else {
#		unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-", fitDistr, "fit-llqqplot.eps", sep = "" ) );
#		warning( "Error while drawing Log-Log Q-Q Plot distribution: '", fitDistr, "'." );
#		warning( "Error message is: ", geterrmessage() );
#		cat( "(ERROR) Log-Log Q-Q Plot failed.", "\n" );
#	}
#	# END Log-Log Q-Q Plot

	# BEGIN P-P Plot
	mg_debug.message( dbg, "[MG::APP::OVERALL::TERAGRID::IAT>> Drawing P-P Plot for '", fitDistr, "' distribution" );
	cat( "## P-P Plot:", "\n" );
	postscript(
		file = paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-", fitDistr, "fit-ppplot.eps", sep = "" ),
		horizontal = FALSE,
		onefile = FALSE
	);
	res <- try(
		do.call(
			"mg_plot_pp",
			c(
				list(iat),
				list(fitDistr),
				fit[[ fitDistr ]]$estimate
			)
		)
	);
	dev.off();
	if ( !is( res, "try-error" ) )
	{
		cat( "\t", paste( "Slope: ", res$slope, " - Intercept: ", res$int ), "\n" );
		cat( "\t", "Correlation Coeff.: ", paste( "Pearson: ", res$cor.pearson$estimate, " @ [", res$cor.pearson$conf.int[1], ",", res$cor.pearson$conf.int[2], "] (p-value: ", res$cor.pearson$p.value, ")", sep = "" ), "\n" );
		cat( "\t", "Area Difference: Absolute: ", res$abs.area.diff, ", Relative: ", res$rel.area.diff, "\n" );
	} else {
		unlink( paste( MG_APP_CONF_TERAGRID_IMGPATH, "/teragrid-overall-iat-", fitDistr, "fit-ppplot.eps", sep = "" ) );
		warning( "Error while drawing P-P Plot distribution: '", fitDistr, "'." );
		warning( "Error message is: ", geterrmessage() );
		cat( "(ERROR) P-P Plot failed.", "\n" );
	}
	# END P-P Plot
}

##@} Performs GoF tests

##@} Distribution fitting (parameters estimated from data set)

sink( file = NULL );

options( warn = oldwarn );
